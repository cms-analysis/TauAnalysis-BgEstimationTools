#include "TauAnalysis/BgEstimationTools/interface/templateBgEstFitAuxFunctions.h"

#include <RooArgList.h>
#include <RooAbsArg.h>
#include <RooRealVar.h>

#include <TRandom3.h>
#include <TMath.h>

#include <iostream>
#include <iomanip>

TRandom3 gRndNum;

double getSampledPull(double pullRMS, double pullMin, double pullMax)
{
  double fluctPull = 0.;
  bool fluctPull_isValid = false;

  while ( !fluctPull_isValid ) {
    double x = gRndNum.Gaus(0., pullRMS);
    if ( x >= pullMin && x <= pullMax ) {
      fluctPull = x;
      fluctPull_isValid = true;
    }
  }

  return fluctPull;
}

void sampleHistogram_stat(const TH1* origHistogram, TH1* fluctHistogram)
{
//--- fluctuate bin-contents by Gaussian distribution
//    with zero mean and standard deviation given by bin-errors

  //std::cout << "<sampleHistogram_stat>:" << std::endl;

  int numBins = origHistogram->GetNbinsX();
  for ( int iBin = 0; iBin < (numBins + 2); ++iBin ) {
    double origBinContent = origHistogram->GetBinContent(iBin);
    double origBinError = origHistogram->GetBinError(iBin);

    double fluctPull = getSampledPull(1., -5., +5.);
    double fluctBinContent = origBinContent + fluctPull*origBinError;

    //std::cout << "iBin = " << iBin << ": origBinContent = " << origBinContent << ","
    //          << " fluctPull = " << fluctPull
    //	        << " --> fluctBinContent = " << fluctBinContent << std::endl;
    
    fluctHistogram->SetBinContent(iBin, fluctBinContent);
    fluctHistogram->SetBinError(iBin, origBinError);
  }
}

void sampleHistogram_sys(TH1* fluctHistogram, const TH1* sysHistogram, 
			 double pullRMS, double pullMin, double pullMax, 
			 int fluctMode)
{
  //std::cout << "<sampleHistogram_sys>:" << std::endl;

  assert(fluctMode == kCoherent || fluctMode == kIncoherent);

//--- fluctuate bin-contents by Gaussian distribution
//    with zero mean and standard deviation given by bin-errors
  double sampledPull = getSampledPull(pullRMS, pullMin, pullMax);

  int numBins = fluctHistogram->GetNbinsX();
  for ( int iBin = 0; iBin < (numBins + 2); ++iBin ) {
    double fluctBinContent = fluctHistogram->GetBinContent(iBin);
    double fluctBinError = fluctHistogram->GetBinError(iBin);
    
    double sysBinContent = sysHistogram->GetBinContent(iBin);
    double sysBinError = sysHistogram->GetBinError(iBin);

    double modBinContent = fluctBinContent + sampledPull*sysBinContent;
    double modBinError = TMath::Sqrt(fluctBinError*fluctBinError + (sampledPull*sysBinError)*(sampledPull*sysBinError));

    //std::cout << "iBin = " << iBin << ": fluctBinContent = " << fluctBinContent << ","
    //	        << " sysBinContent = " << sysBinContent << ", sampledPull = " << sampledPull
    //	        << " --> modBinContent = " << modBinContent << std::endl;

    fluctHistogram->SetBinContent(iBin, modBinContent);
    fluctHistogram->SetBinError(iBin, modBinError);

    if ( fluctMode == kIncoherent ) sampledPull = getSampledPull(pullRMS, pullMin, pullMax);
  }
}

TArrayD getBinning(const TH1* histogram)
{
  TArrayD binning;

//--- check if histogram given as function argument
//    has been created with array of explicitely specified bin-edges
//    or by specifying (numBins, xMin, xMax)
  if ( histogram->GetXaxis()->GetXbins() && histogram->GetXaxis()->GetXbins()->GetSize() > 0 ) {
    binning = (*histogram->GetXaxis()->GetXbins());
  } else {
    int numBins = histogram->GetXaxis()->GetNbins();
    double xMin = histogram->GetXaxis()->GetXmin();
    double xMax = histogram->GetXaxis()->GetXmax();
    double binWidth = (xMax - xMin)/numBins;
    binning.Set(numBins + 1);
    for ( int iBin = 0; iBin < (numBins + 1); ++iBin ) {
      binning[iBin] = xMin + iBin*binWidth;
    }
  }
  
  return binning;
}

double getIntegral(const TH1* histogram)
{
  int numBins = histogram->GetNbinsX();
  return histogram->Integral(0, numBins + 1);
}

double getIntegral(const TH1* histogram, double xMin, double xMax)
{
//--- return sum of entries in histogram given as function argument within range xMin...xMax

  double integral = 0.;

  int numBins = histogram->GetNbinsX();
  for ( int iBin = 0; iBin < (numBins + 2); ++iBin ) {
    double binCenter = histogram->GetBinCenter(iBin);
    
    if ( !(binCenter > xMin && binCenter < xMax) ) continue;

    double binContent = histogram->GetBinContent(iBin);

    integral += binContent;
  }

  return integral;
}


void makeHistogramPositive(TH1* fluctHistogram)
{
  int numBins = fluctHistogram->GetNbinsX();
  for ( int iBin = 0; iBin < (numBins + 2); ++iBin ) {
    if ( fluctHistogram->GetBinContent(iBin) < 0. ) fluctHistogram->SetBinContent(iBin, 0.);
  }
}

TH1* makeConcatenatedHistogram(const std::string& histogramName_concatenated, const std::vector<const TH1*>& histograms, 
			       const std::vector<double_pair>& xRanges)
{
  std::cout << "<makeConcatenatedHistogram>:" << std::endl;

  assert(histograms.size() == xRanges.size());
  
  unsigned numHistograms = histograms.size();

  int numBinsTot = 0;
  
  for ( unsigned iHistogram = 0; iHistogram < numHistograms; ++iHistogram ) {
    const TH1* histogram_i = histograms[iHistogram];
    
    int numBins_i = histogram_i->GetNbinsX();
    for ( int iBin_i = 0; iBin_i < (numBins_i + 2); ++iBin_i ) {
      double binCenter_i = histogram_i->GetBinCenter(iBin_i);
      
      double xMin = xRanges[iHistogram].first;
      double xMax = xRanges[iHistogram].second;

//--- take care that ranges of "original" histograms excluded from fit
//    do not get included in fit of concatenated histogram      
      if ( binCenter_i >= xMin && binCenter_i <= xMax ) ++numBinsTot;
    }
  }

  TH1* histogram_concatenated = new TH1F(histogramName_concatenated.data(), histogramName_concatenated.data(), numBinsTot, -0.5, numBinsTot - 0.5);

  int iBin_concat = 0;

  for ( unsigned iHistogram = 0; iHistogram < numHistograms; ++iHistogram ) {
    const TH1* histogram_i = histograms[iHistogram];

    int numBins_i = histogram_i->GetNbinsX();
    for ( int iBin_i = 0; iBin_i < (numBins_i + 2); ++iBin_i ) {
      double binCenter_i = histogram_i->GetBinCenter(iBin_i);

      double xMin = xRanges[iHistogram].first;
      double xMax = xRanges[iHistogram].second;

//--- take care that ranges of "original" histograms excluded from fit
//    do not get included in fit of concatenated histogram
      if ( binCenter_i >= xMin && binCenter_i <= xMax ) {
	double binContent_i = histogram_i->GetBinContent(iBin_i);
	double binError_i = histogram_i->GetBinError(iBin_i);

	histogram_concatenated->SetBinContent(iBin_concat, binContent_i);
	histogram_concatenated->SetBinError(iBin_concat, binError_i);

	++iBin_concat;
      }
    }
  }

//--- normalize concatenated histogram to total number of entries in first histogram
//    (the same event enters concatenated histograms multiple times,
//     and would hence be "double-counted" if the normalization of the concatenated histogram
//     is not corrected for accordingly)
  double integral_concatenated = getIntegral(histogram_concatenated);
  double norm = getIntegral(histograms[0]);
  if ( integral_concatenated > 0. ) {
    std::cout << "--> scaling histogram_concatenated by factor = " << (norm/integral_concatenated) << std::endl;
    histogram_concatenated->Scale(norm/integral_concatenated);
  }

  return histogram_concatenated;
}

