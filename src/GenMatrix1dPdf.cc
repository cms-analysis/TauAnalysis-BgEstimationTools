#include "TauAnalysis/BgEstimationTools/interface/GenMatrix1dPdf.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

const double epsilon = 1.e-3;

GenMatrix1dPdf::GenMatrix1dPdf(const char* name, const char* title, RooAbsReal& x, 
			       Int_t numBins, Double_t* binBoundaries, TCollection& binProbabilities)
  : RooAbsPdf(name, title),
    x_(std::string(x.GetName()).append("_proxy").data(), "Proxy", this, x)
{
  //std::cout << "<GenMatrix1dPdf::GenMatrix1dPdf>:" << std::endl;
  //std::cout << " name = " << name << std::endl;
  //std::cout << " variable = " << x.GetName() << std::endl;
  //std::cout << " numBins = " << numBins << std::endl;

  if ( numBins != binProbabilities.GetEntries() ) {
    edm::LogError ("GenMatrix1dPdf") 
      << " Number of bin-Probabilities = " << binProbabilities.GetEntries() 
      << " for variable x = " << name << " does not match number of Bins = " << numBins << " !!";
    assert(0);
  }

  binBoundaries_ = new Double_t[numBins + 1];
  for ( Int_t iBin = 0; iBin <= numBins; ++iBin ) {
    //std::cout << " binBoundary[" << iBin << "] = " << binBoundaries[iBin] << std::endl;
    binBoundaries_[iBin] = binBoundaries[iBin];
  }

  binBoundaries_[0] -= epsilon;
  binBoundaries_[numBins] += epsilon;
  
  TIter it(&binProbabilities); 
  while ( RooAbsReal* binProbability = dynamic_cast<RooAbsReal*>(it.Next()) ) {
    RooRealProxy* binProbability_proxy = new RooRealProxy(binProbability->GetName(), "Proxy", this, *binProbability);
    //std::cout << " binProbability[" << binProbabilities_.GetEntries() << "] = " << (Double_t)(*binProbability_proxy) << std::endl;

    if ( (Double_t)(*binProbability_proxy) < 0 ) {
      edm::LogError ("GenMatrix1dPdf") 
	<< " Invalid binProbability[" << binProbabilities_.GetEntries() << ": value must not be negative !!";
      assert(0);
    }

    binProbabilities_.Add(binProbability_proxy);
  }
}

GenMatrix1dPdf::GenMatrix1dPdf(const GenMatrix1dPdf& bluePrint, const char* newName)
  : RooAbsPdf(bluePrint, newName),
    x_(std::string(bluePrint.x_.GetName()).append("_proxy").append("_cloned").data(), this, bluePrint.x_)
{
  Int_t numBins = bluePrint.binProbabilities_.GetEntries();
  binBoundaries_ = new Double_t[numBins + 1];
  for ( Int_t iBin = 0; iBin <= numBins; ++iBin ) {
    binBoundaries_[iBin] = bluePrint.binBoundaries_[iBin];
  }

  TIter it(&bluePrint.binProbabilities_); 
  while ( RooRealProxy* binProbability_proxy = dynamic_cast<RooRealProxy*>(it.Next()) ) {
    binProbabilities_.Add(new RooRealProxy(*binProbability_proxy));
  }
}

GenMatrix1dPdf::~GenMatrix1dPdf()
{
  delete binBoundaries_;

  TIter binProbability(&binProbabilities_); 
  while ( RooRealProxy* it = dynamic_cast<RooRealProxy*>(binProbability.Next()) ) {
    delete it;
  }
}

Double_t GenMatrix1dPdf::evaluate() const
{
  //std::cout << "<GenMatrix1dPdf::evaluate>:" << std::endl;
  //std::cout << " name = " << GetName() << std::endl;
  //std::cout << " variable = " << x_.GetName() << std::endl;
  //std::cout << " x = " << (Double_t)x_ << std::endl;

  Int_t numBins = binProbabilities_.GetEntries();
  for ( Int_t iBin = 0; iBin < numBins; ++iBin ) {
    if ( (Double_t)x_ >= binBoundaries_[iBin] && (Double_t)x_ < binBoundaries_[iBin + 1] ) {
      RooRealProxy* binProbability_proxy = dynamic_cast<RooRealProxy*>(binProbabilities_.At(iBin));
      //std::cout << "--> iBin = " << iBin << ": p = " << (Double_t)(*binProbability_proxy) << std::endl;
      
      return (Double_t)(*binProbability_proxy);
    }
  }

//--- value of x to be evaluated outside 
//    of validity range xMin..xMax of pdf
  edm::LogWarning ("GenMatrix1dPdf::evaluate") 
    << " Value of x = " << x_ << " to be evaluated " 
    << " outside of Range [xMin..xMax[ = [" << binBoundaries_[0] << ".." << binBoundaries_[numBins] << "[" 
    << " of validity of PDF for variable x = " << this->GetName() << " !!";
  return 0.;
}

