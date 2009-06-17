#include "TauAnalysis/BgEstimationTools/interface/DataBinning.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TauAnalysis/BgEstimationTools/interface/binningAuxFunctions.h"

#include <TPRegexp.h>
#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TMath.h>

#include <iostream>
#include <iomanip>
 
DataBinning::DataBinning()
{
  numBins_ = 0;
}

DataBinning::DataBinning(const edm::ParameterSet& cfg) 
  : BinningBase(cfg)
{
  //std::cout << "<DataBinning::DataBinning>:" << std::endl;

  numBins_ = binGrid_->numBins();
  //std::cout << "numBins = " << numBins_ << std::endl;

  binContents_.resize(numBins_);
  binSumw2_.resize(numBins_);
  for ( unsigned iBin = 0; iBin < numBins_; ++iBin ) {
    binContents_[iBin] = 0.;
    binSumw2_[iBin] = 0.;
  }
}

DataBinning::~DataBinning()
{
//--- nothing to be done yet...
}

void DataBinning::bin(const std::vector<double>& x, double weight)
{
  int iBin = binGrid_->binNumber(x);
  if ( iBin >= 0 && iBin < (int)numBins_ ) {
    binContents_[iBin] += weight;
    binSumw2_[iBin] += weight*weight;
  }
}

void DataBinning::print(std::ostream& stream) const
{
  stream << "<DataBinning::print>:" << std::endl;
  stream << " name = " << name_ << std::endl;
  
  const std::vector<std::string>& objVarNames = binGrid_->objVarNames();
  
  double binContent_sum = 0.;
  double binSumw2_sum = 0.;
  for ( unsigned iBin = 0; iBin < numBins_; ++iBin ) {
    stream << " bin " << std::setw(2) << iBin << " (center: ";

    vdouble binCenter = binGrid_->binCenter(iBin);
    if ( binCenter.size() != objVarNames.size() ) {
      edm::LogError ("DataBinning::print") << "Invalid dimension of bin-center vector !!";
      return;
    }

    unsigned numObjVarNames = objVarNames.size();
    for ( unsigned iObjVar = 0; iObjVar < numObjVarNames; ++iObjVar ) {
      stream << objVarNames[iObjVar] << " = " << std::setprecision(3) << std::fixed << binCenter[iObjVar];
      if ( iObjVar < (numObjVarNames - 1) ) stream << ", ";
    }

    stream << "): " << std::setprecision(3) << std::fixed << binContents_[iBin] << " +/- " << TMath::Sqrt(binSumw2_[iBin]) << std::endl;
    
    binContent_sum += binContents_[iBin];
    binSumw2_sum += binSumw2_[iBin];
  }

  stream << " sum: " << std::setprecision(3) << std::fixed << binContent_sum << " +/- " << TMath::Sqrt(binSumw2_sum) << std::endl;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

std::vector<std::string> DataBinning::encodeStringRep() const
{
  std::vector<std::string> buffer = BinningBase::encodeStringRep();

  for ( unsigned iBin = 0; iBin < numBins_; ++iBin ) {
    std::ostringstream meName_binContent;
    meName_binContent << "binContent_region" << (iBin + 1);
    std::ostringstream meValue_binContent;
    meValue_binContent << std::setprecision(3) << std::fixed << binContents_[iBin];
    std::string entry_binContent = encodeBinningStringRep(meName_binContent.str(), "float", meValue_binContent.str());
    buffer.push_back(entry_binContent);
    
    std::ostringstream meName_binSumw2;
    meName_binSumw2 << "binSumw2_region" << (iBin + 1);
    std::ostringstream meValue_binSumw2;
    meValue_binSumw2 << std::setprecision(3) << std::fixed << binSumw2_[iBin];
    std::string entry_binSumw2 = encodeBinningStringRep(meName_binSumw2.str(), "float", meValue_binSumw2.str());
    buffer.push_back(entry_binSumw2);
  }

  return buffer;
}

void DataBinning::decodeStringRep(std::vector<std::string>& buffer)
{
  BinningBase::decodeStringRep(buffer);

  TPRegexp regexpParser_numBins("numBins");
  TPRegexp regexpParser_binContents_entry("binContent_region[[:digit:]]+");
  TPRegexp regexpParser_binContents_binNumber("binContent_region([[:digit:]]+)");
  TPRegexp regexpParser_binSumw2_entry("binSumw2_region[[:digit:]]+");
  TPRegexp regexpParser_binSumw2_binNumber("binSumw2_region([[:digit:]]+)");

  bool numBins_initialized = false;
  std::vector<bool> binContents_initialized;
  std::vector<bool> binSumw2_initialized;
  
  for ( std::vector<std::string>::const_iterator entry = buffer.begin();
	entry != buffer.end(); ++entry ) {
    std::string meName, meType, meValue;
    int error = 0;
    decodeBinningStringRep(*entry, meName, meType, meValue, error);

    if ( error ) {
      edm::LogError ("operator>>") << " Error in parsing string = " << (*entry) << " --> skipping !!";
      continue;
    }

    TString meName_tstring = meName.data();
    
    bool binNumber_error = false;

    if ( regexpParser_numBins.Match(meName_tstring) == 1 ) {
      unsigned numBins = (unsigned)atoi(meValue.data());
      numBins_ = numBins;
      binContents_.resize(numBins_);
      binSumw2_.resize(numBins_);
      binContents_initialized.resize(numBins_);
      binSumw2_initialized.resize(numBins_);
      numBins_initialized = true;
    } else if ( regexpParser_binContents_entry.Match(meName_tstring) == 1 ) {
      if ( !numBins_initialized ) {
	edm::LogError ("operator>>") << " Need to initialize numBins before setting binContents !!";
	continue;
      }
      
      TObjArray* subStrings = regexpParser_binContents_binNumber.MatchS(meName_tstring);
      if ( subStrings->GetEntries() == 2 ) {
	unsigned binNumber = (unsigned)atoi(((TObjString*)subStrings->At(1))->GetString().Data());
	float binContent = atof(meValue.data());
	
	if ( binNumber < numBins_ ) {
	  binContents_[binNumber] = binContent;
	  binContents_initialized[binNumber] = true;
	} else {
	  edm::LogError ("operator>>") << " Bin number = " << binNumber << " decoded from meName = " << meName
				       << " not within numBins = " << numBins_ << " range of binning object !!";
	  continue;
	}
      } else {
	binNumber_error = true;
      }
    } else if ( regexpParser_binSumw2_entry.Match(meName_tstring) == 1 ) {
      if ( !numBins_initialized ) {
	edm::LogError ("operator>>") << " Need to initialize numBins before setting binSumw2 !!";
	continue;
      }

      TObjArray* subStrings = regexpParser_binSumw2_binNumber.MatchS(meName_tstring);
      if ( subStrings->GetEntries() == 2 ) {
	unsigned binNumber = (unsigned)atoi(((TObjString*)subStrings->At(1))->GetString().Data());
	float binSumw2 = atof(meValue.data());

	if ( binNumber < numBins_ ) {
	  binSumw2_[binNumber] = binSumw2;
	  binSumw2_initialized[binNumber] = true;
	} else {
	  edm::LogError ("operator>>") << " Bin number = " << binNumber << " decoded from meName = " << meName
				       << " not within numBins = " << numBins_ << " range of binning object !!";
	  continue;
	}
      } else {
	binNumber_error = true;
      }
    }

    if ( binNumber_error ) {
      edm::LogError ("operator>>") << " Failed to decode bin number from meName = " << meName << " !!";
      continue;
    }
  }

//--- check that all data-members of DataBinning object 
//    have been initialized
  if ( numBins_initialized ) {
    for ( unsigned iBin = 0; iBin < numBins_; ++iBin ) {
      if ( !binContents_initialized[iBin] ) edm::LogError ("operator>>") << " Failed to decode binContents[" << iBin << "] !!";
      if ( !binSumw2_initialized[iBin] ) edm::LogError ("operator>>") << " Failed to decode binSumw2[" << iBin << "] !!";
    }
  } else {
    edm::LogError ("operator>>") << " Failed to decode numBins !!";
  }
}


