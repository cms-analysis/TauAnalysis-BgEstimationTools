#include "TauAnalysis/BgEstimationTools/interface/GenMatrixBgEst1dPdf.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

const float epsilon = 1.e-3;

GenMatrixBgEst1dPdf::GenMatrixBgEst1dPdf(const char* name, const char* title, RooAbsReal& x, 
					 Int_t numBins, Float_t* binBoundaries, TCollection& binProbabilities)
  : RooAbsPdf(name, title),
    x_("x", "Proxy", this, x)
{
  //std::cout << "<GenMatrixBgEst1dPdf::GenMatrixBgEst1dPdf>:" << std::endl;
  //std::cout << " name = " << name << std::endl;
  //std::cout << " variable = " << x.GetName() << std::endl;
  //std::cout << " numBins = " << numBins << std::endl;

  if ( numBins != binProbabilities.GetEntries() ) {
    edm::LogError ("GenMatrixBgEst1dPdf") << " Number of bin-Probabilities = " << binProbabilities.GetEntries() 
					  << " for variable x = " << name << " does not match number of Bins = " << numBins << " !!";
    assert(0);
  }

  binBoundaries_ = new Float_t[numBins + 1];
  for ( Int_t iBin = 0; iBin <= numBins; ++iBin ) {
    std::cout << " binBoundary[" << iBin << "] = " << binBoundaries[iBin] << std::endl;
    binBoundaries_[iBin] = binBoundaries[iBin];
  }

  binBoundaries_[0] -= epsilon;
  binBoundaries_[numBins] += epsilon;
  
  TIter it(&binProbabilities); 
  while ( RooAbsReal* binProbability = dynamic_cast<RooAbsReal*>(it.Next()) ) {
    RooRealProxy* binProbability_proxy = new RooRealProxy(binProbability->GetName(), "Proxy", this, *binProbability);
    binProbabilities_.Add(binProbability_proxy);
  }
}

GenMatrixBgEst1dPdf::GenMatrixBgEst1dPdf(const GenMatrixBgEst1dPdf& bluePrint, const char* newName)
  : RooAbsPdf(bluePrint, newName),
    x_("x", this, bluePrint.x_)
{
  Int_t numBins = bluePrint.binProbabilities_.GetEntries();
  binBoundaries_ = new Float_t[numBins + 1];
  for ( Int_t iBin = 0; iBin <= numBins; ++iBin ) {
    binBoundaries_[iBin] = bluePrint.binBoundaries_[iBin];
  }

  TIter it(&bluePrint.binProbabilities_); 
  while ( RooRealProxy* binProbability_proxy = dynamic_cast<RooRealProxy*>(it.Next()) ) {
    binProbabilities_.Add(new RooRealProxy(*binProbability_proxy));
  }
}

GenMatrixBgEst1dPdf::~GenMatrixBgEst1dPdf()
{
  delete binBoundaries_;

  TIter binProbability(&binProbabilities_); 
  while ( RooRealProxy* it = dynamic_cast<RooRealProxy*>(binProbability.Next()) ) {
    delete it;
  }
}

Double_t GenMatrixBgEst1dPdf::evaluate() const
{
  Int_t numBins = binProbabilities_.GetEntries();
  for ( Int_t iBin = 0; iBin < numBins; ++iBin ) {
    if ( x_ >= binBoundaries_[iBin] && x_ < binBoundaries_[iBin + 1] ) {
      RooRealProxy* binProbability_proxy = dynamic_cast<RooRealProxy*>(binProbabilities_.At(iBin));

      return *binProbability_proxy;
    }
  }

//--- value of x to be evaluated outside 
//    of validity range xMin..xMax of pdf
  edm::LogWarning ("GenMatrixBgEst1dPdf::evaluate") << " Value of x = " << x_ << " to be evaluated " 
						    << " outside of Range [xMin..xMax[ = [" 
						    << binBoundaries_[0] << ".." << binBoundaries_[numBins] << "[" 
						    << " of validity of PDF for variable x = " << this->GetName() << " !!";
  return 0.;
}

