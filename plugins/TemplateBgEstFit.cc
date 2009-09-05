#include "TauAnalysis/BgEstimationTools/plugins/TemplateBgEstFit.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "TauAnalysis/DQMTools/interface/dqmAuxFunctions.h"
#include "TauAnalysis/DQMTools/interface/generalAuxFunctions.h"

#include "TauAnalysis/BgEstimationTools/interface/histogramAuxFunctions.h"
#include "TauAnalysis/BgEstimationTools/interface/BgEstMean.h"
#include "TauAnalysis/BgEstimationTools/interface/BgEstMedian.h"
#include "TauAnalysis/BgEstimationTools/interface/BgEstCovMatrix.h"

#include <TCanvas.h>
#include <TROOT.h>
#include <TMath.h>
#include <TEllipse.h>
#include <TColor.h>
#include <TMarker.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TArrayD.h>
#include <TLegend.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TStyle.h>

#include <RooAddPdf.h>
#include <RooGaussian.h>
#include <RooProdPdf.h>
#include <RooArgSet.h>
#include <RooGlobalFunc.h>
#include <RooPlot.h>
#include <RooFit.h>
#include <RooLinkedList.h>
#include <RooCmdArg.h>
#include <RooTFnPdfBinding.h>

#include <iostream>

const int defaultCanvasSizeX = 800;
const int defaultCanvasSizeY = 600;

const double maxNorm = 1.e+6;
const double normStartValue = 1.e+3;

const std::string fluctMode_coherent = "coherent";
const std::string fluctMode_incoherent = "incoherent";

enum { kCoherent, kIncoherent };

const std::string fitMode_1dPlus = "1d+";
const std::string fitMode_Nd = "Nd";

enum { k1dPlus, kNd };

TRandom3 gRndNum;

const double epsilon = 1.e-3;
const double epsilon_integral = 5.e-2;

const int fitStatus_converged = 0;

typedef std::pair<double, double> double_pair;

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
//--- return total number of entries in histogram given as function argument
//    (including underflow and overflow bins)

  if ( histogram->InheritsFrom(TH3::Class()) ) {
    const TH3* histogram3d = dynamic_cast<const TH3*>(histogram);
    assert(histogram3d);
    int numBinsX = histogram3d->GetNbinsX();
    int numBinsY = histogram3d->GetNbinsY();
    int numBinsZ = histogram3d->GetNbinsZ();
    return histogram3d->Integral(0, numBinsX + 1, 0, numBinsY + 1, 0, numBinsZ + 1);
  } else if ( histogram->InheritsFrom(TH2::Class()) ) {
    const TH2* histogram2d = dynamic_cast<const TH2*>(histogram);
    assert(histogram2d);
    int numBinsX = histogram2d->GetNbinsX();
    int numBinsY = histogram2d->GetNbinsY();
    return histogram2d->Integral(0, numBinsX + 1, 0, numBinsY + 1);
  } else {
    int numBins = histogram->GetNbinsX();
    return histogram->Integral(0, numBins + 1);
  }
}

double getIntegral1d(const TH1* histogram, double xMin, double xMax)
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

TH1* makeConcatenatedHistogram(const std::string& concatHistogramName, const std::vector<TH1*>& histograms, 
			       const std::vector<double_pair>& xRanges)
{
  std::cout << "<makeConcatenatedHistogram>:" << std::endl;

  unsigned numHistograms = histograms.size();

  int numBinsTot = 0;

  for ( unsigned iHistogram = 0; iHistogram < numHistograms; ++iHistogram ) {
    const TH1* histogram_i = histograms[iHistogram];
    numBinsTot += (histogram_i->GetNbinsX() + 2);
  }

  TH1* concatHistogram = new TH1F(concatHistogramName.data(), concatHistogramName.data(), numBinsTot, -0.5, numBinsTot - 0.5);

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

	concatHistogram->SetBinContent(iBin_concat, binContent_i);
	concatHistogram->SetBinError(iBin_concat, binError_i);
      } else {
	concatHistogram->SetBinContent(iBin_concat, 0.);
	concatHistogram->SetBinError(iBin_concat, 0.);
      }

      ++iBin_concat;
    }
  }

//--- normalize concatenated histogram to total number of entries in first histogram
//    (the same event enters concatenated histograms multiple times,
//     and would hence be "double-counted" if the normalization of the concatenated histogram
//     is not corrected for accordingly)
  double integral_concatenated = getIntegral(concatHistogram);
  double norm = getIntegral(histograms[0]);
  if ( integral_concatenated > 0. ) {
    std::cout << "--> scaling concatHistogram by factor = " << (norm/integral_concatenated) << std::endl;
    concatHistogram->Scale(norm/integral_concatenated);
  }

  return concatHistogram;
}

void unpackFitResult(const RooFitResult* fitResult, TVectorD& mean, TMatrixD& cov)
{
  const RooArgList& fitParameter = fitResult->floatParsFinal();
  int numFitParameter = fitParameter.getSize();
  mean.ResizeTo(numFitParameter);
  cov.ResizeTo(numFitParameter, numFitParameter);
  for ( int iX = 0; iX < numFitParameter; ++iX ) {
    const RooAbsArg* paramX_arg = fitParameter.at(iX);
    const RooRealVar* paramX = dynamic_cast<const RooRealVar*>(paramX_arg);
    assert(paramX != NULL);
    mean(iX) = paramX->getVal();
    double sigmaX = paramX->getError();    
    for ( int iY = 0; iY < numFitParameter; ++iY ) {
      if ( iY == iX ) {
	cov(iX, iX) = sigmaX*sigmaX;
      } else {
	const RooAbsArg* paramY_arg = fitParameter.at(iY);
	const RooRealVar* paramY = dynamic_cast<const RooRealVar*>(paramY_arg);
	assert(paramY != NULL);
	double sigmaY = paramY->getError();
	double corrXY = fitResult->correlation(*paramX_arg, *paramY_arg);
	cov(iX, iY) = sigmaX*sigmaY*corrXY;
      }
    }
  }
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

TemplateBgEstFit::dataDistr1dType::dataDistr1dType(const std::string& processName, const std::string& varName, 
						   const std::string& meName, RooRealVar* x, bool cutUnfittedRegion)
  : processName_(processName),
    varName_(varName),
    meName_(meName),
    xRef_(x),
    cutUnfittedRegion_(cutUnfittedRegion),
    histogram_(0),
    dataHist_(0), 
    error_(0)
{}

TemplateBgEstFit::dataDistr1dType::~dataDistr1dType()
{
  delete histogram_;
  delete dataHist_;
  delete fluctHistogram_;
}

void TemplateBgEstFit::dataDistr1dType::initialize()
{
  //std::cout << "<dataDistr1dType::initialize>:" << std::endl;

  DQMStore& dqmStore = (*edm::Service<DQMStore>());
  me_ = dqmStore.get(meName_);
  if ( !me_ ) {
    edm::LogError ("dataDistr1dType") << " Failed to access dqmMonitorElement = " << meName_ << " !!";
    error_ = 1;
    return;
  }

//--- make sure that x-axis range excluded from fit
//    does really have no effect on fit results
//    (differences between the event selection criteria
//     applied in the final analysis and the background enriched samples
//     may cause sizeable deviations between the templates obtained from the background enriched samples 
//     and the distributions observed in the final analysis;
//     even in case regions with deviations are not within the fitted region,
//     they may affect the fit results via differences in the fractions of events 
//     that are outside of the region included in the fit,
//     because the normalization of PDFs depends on these fractions)
  histogram_ = (TH1*)me_->getTH1()->Clone();
  if ( cutUnfittedRegion_ ) {
    const RooAbsBinning& xRange = xRef_->getBinning();
    double xMin = xRange.lowBound();
    double xMax = xRange.highBound();
    int numBins = histogram_->GetNbinsX();
    for ( int iBin = 0; iBin < (numBins + 2); ++iBin ) {
      double binCenter = histogram_->GetBinCenter(iBin);
      
      if ( !(binCenter > xMin && binCenter < xMax) ) {
	histogram_->SetBinContent(iBin, 0.);
	histogram_->SetBinError(iBin, 0.);
      }
    }
  }

  fluctHistogram_ = (TH1*)histogram_->Clone();

  dataHistName_ = std::string(processName_).append("_").append(varName_).append("_rooDataHist");
  buildFitData();
}

void TemplateBgEstFit::dataDistr1dType::buildFitData()
{
  //std::cout << "<dataDistr1dType::buildFitData>:" << std::endl;
  //std::cout << " dataHistName = " << dataHistName_ << std::endl;
  dataHist_ = new RooDataHist(dataHistName_.data(), dataHistName_.data(), *xRef_, fluctHistogram_);
  //std::cout << " dataHist = " << dataHist_ << std::endl;
}

void TemplateBgEstFit::dataDistr1dType::fluctuate(bool, bool)
{  
  sampleHistogram_stat(histogram_, fluctHistogram_);
  
  makeHistogramPositive(fluctHistogram_);
  
  //std::cout << "<dataDistr1dType::fluctuate>:" << std::endl;
  //std::cout << "--> deleting dataHist..." << std::endl;
  //std::cout << " dataHist = " << dataHist_ << std::endl;

  delete dataHist_;
  buildFitData();
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

TemplateBgEstFit::dataDistrNdType::dataDistrNdType(int fitMode, bool cutUnfittedRegion)
  : numDimensions_(0),
    fitMode_(fitMode),
    cutUnfittedRegion_(cutUnfittedRegion),
    fitData_(0),
    normCorrFactor_(1.),
    auxHistogram_(0),
    error_(0)
{}

TemplateBgEstFit::dataDistrNdType::~dataDistrNdType()
{
  //std::cout << "<dataDistrNdType::~dataDistrNdType>:" << std::endl;

  for ( std::map<std::string, dataDistr1dType*>::iterator it = dataEntries1d_.begin();
	it != dataEntries1d_.end(); ++it ) {
    delete it->second;
  }

  delete fitData_;
  
  delete auxHistogram_;
}

void TemplateBgEstFit::dataDistrNdType::addElement(const std::string& varName, RooRealVar* x, const std::string& meName)
{
  varNames_.push_back(varName);
  
  dataDistr1dType* dataEntry1d = new dataDistr1dType("data", varName, meName, x, cutUnfittedRegion_);
  if ( dataEntry1d->error_ ) error_ = 1;
  dataEntries1d_[varName] = dataEntry1d;
  
  ++numDimensions_;
}

void TemplateBgEstFit::dataDistrNdType::initialize()
{
  std::cout << "<dataDistrNdType::initialize>:" << std::endl;

  for ( std::map<std::string, dataDistr1dType*>::iterator dataEntry1d = dataEntries1d_.begin();
	dataEntry1d != dataEntries1d_.end(); ++dataEntry1d ) { 
    dataEntry1d->second->initialize();
    if ( dataEntry1d->second->error_ ) error_ = 1;

    std::cout << "processName = " << dataEntry1d->second->processName_ << std::endl;

    const RooAbsBinning& xRange = dataEntry1d->second->xRef_->getBinning();
    double xMin = xRange.lowBound();
    double xMax = xRange.highBound();
    double integral_fitted = getIntegral1d(dataEntry1d->second->histogram_, xMin, xMax);
    std::cout << " integral_fitted = " << integral_fitted << std::endl;
    
    double integral = getIntegral(dataEntry1d->second->me_->getTH1());
    std::cout << " integral = " << integral << std::endl;
    
    if ( integral_fitted > 0. ) normCorrFactor_ *= (integral/integral_fitted);
  }

  std::cout << "--> using normCorrFactor = " << normCorrFactor_ << std::endl;

  buildFitData();
}

void TemplateBgEstFit::dataDistrNdType::buildFitData()
{
  //std::cout << "<dataDistrNdType::buildFitData>:" << std::endl;
  //std::cout << " numDimensions = " << numDimensions_ << std::endl;

  assert(fitMode_ == k1dPlus || fitMode_ == kNd);

  if ( numDimensions_ == 1 ) {
    std::string varName_1 = varNames_.front();
    dataDistr1dType* dataEntry1d = dataEntries1d_[varName_1];

    std::string fitDataName = "fitData_rooDataHist";
    //std::cout << " fitDataName = " << fitDataName << std::endl;

    fitData_ = new RooDataHist(fitDataName.data(), fitDataName.data(), *dataEntry1d->xRef_, dataEntry1d->fluctHistogram_);
    //std::cout << " fitData = " << fitData_ << std::endl;
  } else {
    if ( fitMode_ == k1dPlus ) {
      std::vector<TH1*> histograms;
      std::vector<double_pair> xRanges;
      for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
	std::string varName = varNames_[iDimension];
	dataDistr1dType* dataEntry1d = dataEntries1d_[varName];
	histograms.push_back(dataEntry1d->histogram_);
	const RooAbsBinning& xRange = dataEntry1d->xRef_->getBinning();
	double xMin = xRange.lowBound();
	double xMax = xRange.highBound();
	xRanges.push_back(double_pair(xMin, xMax));
      }

      std::string auxHistogramName = "auxHistogram_1dPlus";
      delete auxHistogram_;
      auxHistogram_ = makeConcatenatedHistogram(auxHistogramName, histograms, xRanges);
    } else if ( fitMode_ == kNd ) {
      if ( numDimensions_ == 2 ) {
	const std::string& varName_1 = varNames_[0];
	TH1* histogram_1 = dataEntries1d_[varName_1]->histogram_;
	int numBins_1 = histogram_1->GetXaxis()->GetNbins();
	TArrayD binning_1 = getBinning(histogram_1);
	double integral_1 = getIntegral(histogram_1);
	  
	const std::string& varName_2 = varNames_[1];
	TH1* histogram_2 = dataEntries1d_[varName_2]->histogram_;
	int numBins_2 = histogram_2->GetXaxis()->GetNbins();
	TArrayD binning_2 = getBinning(histogram_2);
	double integral_2 = getIntegral(histogram_2);

	if ( TMath::Abs(integral_1 - integral_2) > epsilon_integral*0.5*(integral_1 + integral_2) ) {
	  static bool isFirstWarning = true;
	  if ( isFirstWarning )	
	    edm::LogWarning ("buildFitData") << " Difference in integrals between first (" << varName_1 << " = " << integral_1 << ")"
					     << " and second (" << varName_2 << " = " << integral_2 << ") Histogram detected !!";
	  isFirstWarning = false;
	}
	  
	std::string auxHistogramName = std::string(histogram_1->GetName()).append("_").append(histogram_2->GetName());
	delete auxHistogram_;
	auxHistogram_ = new TH2F(auxHistogramName.data(), auxHistogramName.data(), 
				 numBins_1, binning_1.GetArray(), numBins_2, binning_2.GetArray());
	for ( int iBin_1 = 0; iBin_1 < (numBins_1 + 2); ++iBin_1 ) {
	  double binCenter_1 = histogram_1->GetBinCenter(iBin_1);
	  double binContent_1 = histogram_1->GetBinContent(iBin_1)/integral_1;
	  double binError_1 = histogram_1->GetBinError(iBin_1)/integral_1;
	    
	  for ( int iBin_2 = 0; iBin_2 < (numBins_2 + 2); ++iBin_2 ) {
	    double binCenter_2 = histogram_2->GetBinCenter(iBin_2);
	    double binContent_2 = histogram_2->GetBinContent(iBin_2)/integral_2;
	    double binError_2 = histogram_2->GetBinError(iBin_2)/integral_2;
	    
	    double binContent_combined = binContent_1*binContent_2;
	    double binError_d1_2 = binError_1*binContent_2;
	    double binError_1_d2 = binContent_1*binError_2;
	    double binError2_combined = binError_d1_2*binError_d1_2 + binError_1_d2*binError_1_d2;
	    
	    int iBin_aux = auxHistogram_->FindBin(binCenter_1, binCenter_2);
	    //std::cout << " iBin_aux = " << iBin_aux << std::endl;
	    
	    auxHistogram_->SetBinContent(iBin_aux, integral_1*binContent_combined);
	    auxHistogram_->SetBinError(iBin_aux, integral_1*TMath::Sqrt(binError2_combined));

	    //std::cout << "setting binContent(" << varName_1 << " = " << binCenter_1 << ", " << varName_2 << " = " << binCenter_2 << ")"
	    //  	<< " = " << auxHistogram_->GetBinContent(iBin_aux)
	    //          << " +/- binError = " << auxHistogram_->GetBinError(iBin_aux) << std::endl;
	  }
	}

	//std::cout << " integral(auxHistogram) = " << getIntegral(auxHistogram_) << ", integral_1 = " << integral_1 << std::endl;
	assert(TMath::Abs(getIntegral(auxHistogram_) - integral_1) < epsilon_integral*0.5*(getIntegral(auxHistogram_) + integral_1));
      } else if ( numDimensions_ == 3 ) {
	const std::string& varName_1 = varNames_[0];
	TH1* histogram_1 = dataEntries1d_[varName_1]->histogram_;
	int numBins_1 = histogram_1->GetXaxis()->GetNbins();
	TArrayD binning_1 = getBinning(histogram_1);
	double integral_1 = getIntegral(histogram_1);
      
	const std::string& varName_2 = varNames_[1];
	TH1* histogram_2 = dataEntries1d_[varName_2]->histogram_;
	int numBins_2 = histogram_2->GetXaxis()->GetNbins();
	TArrayD binning_2 = getBinning(histogram_2);
	double integral_2 = getIntegral(histogram_2);
      
	const std::string& varName_3 = varNames_[2];
	TH1* histogram_3 = dataEntries1d_[varName_3]->histogram_;
	int numBins_3 = histogram_3->GetXaxis()->GetNbins();
	TArrayD binning_3 = getBinning(histogram_3);
	double integral_3 = getIntegral(histogram_3);
	
	if ( TMath::Abs(integral_1 - integral_2) > epsilon_integral*0.5*(integral_1 + integral_2) ||
	     TMath::Abs(integral_1 - integral_3) > epsilon_integral*0.5*(integral_1 + integral_3) ||
	     TMath::Abs(integral_2 - integral_3) > epsilon_integral*0.5*(integral_2 + integral_3) ) {
	  static bool isFirstWarning = true;
	  if ( isFirstWarning )	
	    edm::LogWarning ("buildFitData") << " Difference in integrals between first (" << varName_1 << " = " << integral_1 << "),"
					     << " second (" << varName_2 << " = " << integral_2 << ")"
					     << " and third (" << varName_3 << " = " << integral_3 << ") Histogram detected !!";
	  isFirstWarning = false;
	}

	std::string auxHistogramName
	  = std::string(histogram_1->GetName()).append("_").append(histogram_2->GetName()).append("_").append(histogram_3->GetName());
	delete auxHistogram_;
	auxHistogram_ = new TH3F(auxHistogramName.data(), auxHistogramName.data(), 
				 numBins_1, binning_1.GetArray(), numBins_2, binning_2.GetArray(), numBins_3, binning_3.GetArray());
	for ( int iBin_1 = 0; iBin_1 < (numBins_1 + 2); ++iBin_1 ) {
	  double binCenter_1 = histogram_1->GetBinCenter(iBin_1);
	  double binContent_1 = histogram_1->GetBinContent(iBin_1)/integral_1;
	  double binError_1 = histogram_1->GetBinError(iBin_1)/integral_1;
	
	  for ( int iBin_2 = 0; iBin_2 < (numBins_2 + 2); ++iBin_2 ) {
	    double binCenter_2 = histogram_2->GetBinCenter(iBin_2);
	    double binContent_2 = histogram_2->GetBinContent(iBin_2)/integral_2;
	    double binError_2 = histogram_2->GetBinError(iBin_2)/integral_2;
	    
	    for ( int iBin_3 = 0; iBin_3 < (numBins_3 + 2); ++iBin_3 ) {
	      double binCenter_3 = histogram_3->GetBinCenter(iBin_3);
	      double binContent_3 = histogram_3->GetBinContent(iBin_3)/integral_3;
	      double binError_3 = histogram_3->GetBinError(iBin_3)/integral_3;
	      
	      double binContent_combined = binContent_1*binContent_2*binContent_3;
	      double binError_d1_2_3 = binError_1*binContent_2*binContent_3;
	      double binError_1_d2_3 = binContent_1*binError_2*binContent_3;
	      double binError_1_2_d3 = binContent_1*binContent_2*binError_3;
	      double binError2_combined 
		= binError_d1_2_3*binError_d1_2_3 + binError_1_d2_3*binError_1_d2_3 + binError_1_2_d3*binError_1_2_d3;
	      
	      int iBin_aux = auxHistogram_->FindBin(binCenter_1, binCenter_2, binCenter_3);
	    
	      auxHistogram_->SetBinContent(iBin_aux, integral_1*binContent_combined);
	      auxHistogram_->SetBinError(iBin_aux, integral_1*TMath::Sqrt(binError2_combined));
	    }
	  }
	}
	
	//std::cout << " integral(auxHistogram) = " << getIntegral(auxHistogram_) << ", integral_1 = " << integral_1 << std::endl;
        assert(TMath::Abs(getIntegral(auxHistogram_) - integral_1) < epsilon_integral*0.5*(getIntegral(auxHistogram_) + integral_1));
      } else {
	edm::LogError ("buildFitData") << " Number of Dimensions = " << numDimensions_ << " exceed maximum" 
				       << " supported for fitMode = " << fitMode_ << " !!";
	error_ = 1;
	return;
      }
    } 

    std::string fitDataName = "fitData_rooDataHist";
    
    TObjArray fitData_varsCollection;
    for ( std::map<std::string, dataDistr1dType*>::iterator dataEntry1d = dataEntries1d_.begin();
	  dataEntry1d != dataEntries1d_.end(); ++dataEntry1d ) {
      fitData_varsCollection.Add(dataEntry1d->second->xRef_);
    }
    
    std::string fitData_varsArgName = std::string("fitData").append("_varsArgs");
    RooArgList fitData_varsArgs(fitData_varsCollection, fitData_varsArgName.data());
    
    fitData_ = new RooDataHist(fitDataName.data(), fitDataName.data(), fitData_varsArgs, auxHistogram_);    
  } 
}

void TemplateBgEstFit::dataDistrNdType::fluctuate(bool, bool)
{
  for ( std::map<std::string, dataDistr1dType*>::iterator dataEntry1d = dataEntries1d_.begin();
	dataEntry1d != dataEntries1d_.end(); ++dataEntry1d ) {
    dataEntry1d->second->fluctuate(true, false);
  }

  //std::cout << "<dataDistrNdType::fluctuate>:" << std::endl;
  //std::cout << "--> deleting fitData..." << std::endl;
  //std::cout << " fitData = " << fitData_ << std::endl;
  delete fitData_;
  buildFitData();
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

TemplateBgEstFit::modelTemplate1dType::modelTemplate1dType(const std::string& processName, const std::string& varName,
							   const std::string& meName, RooRealVar* x, bool cutUnfittedRegion,
							   bool applySmoothing, const edm::ParameterSet& cfgSmoothing)
  : dataDistr1dType(processName, varName, meName, x, cutUnfittedRegion),
    pdf1d_(0),
    applySmoothing_(applySmoothing), cfgSmoothing_(cfgSmoothing),
    auxTF1Wrapper_(0)
{}

TemplateBgEstFit::modelTemplate1dType::~modelTemplate1dType()
{
  delete pdf1d_;
  delete auxTF1Wrapper_;
}

void TemplateBgEstFit::modelTemplate1dType::initialize()
{
  //std::cout << "<modelTemplate1dType::initialize>:" << std::endl;

  dataDistr1dType::initialize();

  if ( error_ ) return;
 
  pdf1dName_ = std::string(processName_).append("_").append(varName_).append("_").append("pdf");
  buildPdf();

  DQMStore& dqmStore = (*edm::Service<DQMStore>());
  for ( std::vector<sysFluctDefType>::iterator sysErrFluctuation = sysErrFluctuations_.begin();
	sysErrFluctuation != sysErrFluctuations_.end(); ++sysErrFluctuation ) {
    sysErrFluctuation->me_ = dqmStore.get(sysErrFluctuation->meName_);
    if ( !sysErrFluctuation->me_ ) {
      edm::LogError ("modelTemplate1dType") << " Failed to access dqmMonitorElement = " << sysErrFluctuation->meName_ << " !!";
      error_ = 1;
      continue;
    }

//--- check that histograms representing systematic uncertainties have the same binning
//    as that representing expectation
    if ( !isCompatibleBinning(histogram_, sysErrFluctuation->me_->getTH1()) ) {
      edm::LogError ("modelTemplate1dType") << " Incompatible binning of histograms " << meName_ 
					    << " and " << sysErrFluctuation->meName_ << " !!";
      error_ = 1;
      continue;
    }
  }
}

void TemplateBgEstFit::modelTemplate1dType::buildPdf()
{
  if ( !applySmoothing_ ) {
    pdf1d_ = new RooHistPdf(pdf1dName_.data(), pdf1dName_.data(), *xRef_, *dataHist_);
  } else {

    bool isFirstFit = (!auxTF1Wrapper_);

    if ( isFirstFit ) {
      std::string pluginTypeTF1Wrapper = cfgSmoothing_.getParameter<std::string>("pluginType");
      delete auxTF1Wrapper_;
      auxTF1Wrapper_ = TF1WrapperPluginFactory::get()->create(pluginTypeTF1Wrapper, cfgSmoothing_);
    } else {
      auxTF1Wrapper_->reinitializeTF1Parameter();
    }

    std::string fitOption = ( isFirstFit ) ? "RB0" : "RB0Q";
    
    fluctHistogram_->Fit(auxTF1Wrapper_->getTF1(), fitOption.data());

    // CMSSW_3_1_x only (because needs ROOT version 5.22)
    pdf1d_ = new RooTFnPdfBinding(pdf1dName_.data(), pdf1dName_.data(), auxTF1Wrapper_->getTF1(), RooArgList(*xRef_));
  } 
}

void TemplateBgEstFit::modelTemplate1dType::fluctuate(bool fluctStat, bool fluctSys)
{
  if ( fluctStat ) {
    sampleHistogram_stat(histogram_, fluctHistogram_);
  } else {
    int numBins = histogram_->GetNbinsX();
    for ( int iBin = 0; iBin < (numBins + 2); ++iBin ) {
      double binContent = histogram_->GetBinContent(iBin);
      double binError = histogram_->GetBinError(iBin);

      fluctHistogram_->SetBinContent(iBin, binContent);
      fluctHistogram_->SetBinError(iBin, binError);
    }
  }

  if ( fluctSys ) {
    for ( std::vector<sysFluctDefType>::const_iterator sysErrFluctuation = sysErrFluctuations_.begin();
	  sysErrFluctuation != sysErrFluctuations_.end(); ++sysErrFluctuation ) {
      sampleHistogram_sys(fluctHistogram_, sysErrFluctuation->me_->getTH1(),  
			  sysErrFluctuation->pullRMS_, sysErrFluctuation->pullMin_, sysErrFluctuation->pullMax_, 
			  sysErrFluctuation->fluctMode_);
    }
  }

  makeHistogramPositive(fluctHistogram_);

  if ( !applySmoothing_ ) {
    //std::cout << "<modelTemplate1dType::fluctuate>:" << std::endl;
    //std::cout << "--> deleting dataHist..." << std::endl;
    //std::cout << " dataHist = " << dataHist_ << std::endl;
    delete dataHist_;
    buildFitData();
  }
  
  delete pdf1d_;
  buildPdf();
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

TemplateBgEstFit::modelTemplateNdType::modelTemplateNdType(const std::string& processName, const edm::ParameterSet& cfgProcess,
							   int fitMode, bool cutUnfittedRegion, bool applyNormConstraint, 
							   double valueNormConstraint, double sigmaNormConstraint)
  : processName_(processName),
    applyNormConstraint_(applyNormConstraint),
    pdfNormConstraint_(0), 
    meanNormConstraint_(0),
    sigmaNormConstraint_(0),
    pdf_(0),
    numDimensions_(0),
    fitMode_(fitMode),
    cutUnfittedRegion_(cutUnfittedRegion),
    normCorrFactor_(1.),
    auxHistogram_(0),
    auxDataHist_(0),
    error_(0)
{
  edm::ParameterSet cfgDrawOptions_process = cfgProcess.getParameter<edm::ParameterSet>("drawOptions");
  lineColor_ = cfgDrawOptions_process.getParameter<int>("lineColor");
  lineStyle_ = cfgDrawOptions_process.getParameter<int>("lineStyle");
  lineWidth_ = cfgDrawOptions_process.getParameter<int>("lineWidth");

  std::string normName = std::string(processName_).append("_").append("norm");
  double norm_initial = ( cfgProcess.exists("norm") ) ? 
    cfgProcess.getParameter<edm::ParameterSet>("norm").getParameter<double>("initial") : normStartValue;
  norm_ = new RooRealVar(normName.data(), normName.data(), norm_initial, 0., maxNorm);

  if ( applyNormConstraint_ ) {
    std::cout << "<modelTemplateNdType>: constraining norm = " << valueNormConstraint << " +/- " << sigmaNormConstraint << ","
	      << " process = " << processName_ << std::endl;
    
    std::string meanNormConstraintName = std::string(processName_).append("_").append("meanNormConstraint");
    meanNormConstraint_ = new RooConstVar(meanNormConstraintName.data(), meanNormConstraintName.data(), valueNormConstraint);
    std::string sigmaNormConstraintName = std::string(processName_).append("_").append("sigmaNormConstraint");
    sigmaNormConstraint_ = new RooConstVar(sigmaNormConstraintName.data(), sigmaNormConstraintName.data(), sigmaNormConstraint);
    
    std::string pdfNormConstraintName = std::string(processName_).append("_").append("pdfNormConstraint");
    pdfNormConstraint_ = new RooGaussian(pdfNormConstraintName.data(), pdfNormConstraintName.data(), 
					 *norm_, *meanNormConstraint_, *sigmaNormConstraint_);
  }
}

TemplateBgEstFit::modelTemplateNdType::~modelTemplateNdType()
{
  for ( std::map<std::string, modelTemplate1dType*>::iterator it = processEntries1d_.begin();
	it != processEntries1d_.end(); ++it ) {
    delete it->second;
  }

  delete auxHistogram_;
  delete auxDataHist_;

  if ( pdfIsOwned_ ) delete pdf_;

  delete norm_;

  delete pdfNormConstraint_;
  delete meanNormConstraint_;
  delete sigmaNormConstraint_;
}

void TemplateBgEstFit::modelTemplateNdType::addElement(const std::string& varName, RooRealVar* x, const std::string& meName,
						       bool applySmoothing, const edm::ParameterSet& cfgSmoothing)
{
  varNames_.push_back(varName);
  
  modelTemplate1dType* processEntry1d = new modelTemplate1dType(processName_, varName, meName, x, cutUnfittedRegion_,
								applySmoothing, cfgSmoothing);
  if ( processEntry1d->error_ ) error_ = 1;
  processEntries1d_[varName] = processEntry1d;
  
  ++numDimensions_;
}

void TemplateBgEstFit::modelTemplateNdType::initialize()
{
  std::cout << "<modelTemplateNdType::initialize>:" << std::endl;

  for ( std::map<std::string, modelTemplate1dType*>::iterator processEntry1d = processEntries1d_.begin();
	processEntry1d != processEntries1d_.end(); ++processEntry1d ) {
    processEntry1d->second->initialize();
    if ( processEntry1d->second->error_ ) error_ = 1;

    std::cout << "processName = " << processEntry1d->second->processName_ << std::endl;
      
    const RooAbsBinning& xRange = processEntry1d->second->xRef_->getBinning();
    double xMin = xRange.lowBound();
    double xMax = xRange.highBound();
    double integral_fitted = getIntegral1d(processEntry1d->second->histogram_, xMin, xMax);
    std::cout << " integral_fitted = " << integral_fitted << std::endl;
    
    double integral = getIntegral(processEntry1d->second->me_->getTH1());
    std::cout << " integral = " << integral << std::endl;
    
    if ( integral_fitted > 0. ) normCorrFactor_ *= (integral/integral_fitted);
  }

  std::cout << "--> using normCorrFactor = " << normCorrFactor_ << std::endl;

  buildPdf();
}

void TemplateBgEstFit::modelTemplateNdType::buildPdf()
{
  assert(fitMode_ == k1dPlus || fitMode_ == kNd);

  if ( numDimensions_ == 1 ) {
    const std::string varName = varNames_.front();
    pdfName_ = processEntries1d_[varName]->pdf1dName_;
    pdf_ = processEntries1d_[varName]->pdf1d_;
    pdfIsOwned_ = false;
  } else {
    if ( fitMode_ == k1dPlus ) {
      std::vector<TH1*> histograms;
      std::vector<double_pair> xRanges;
      for ( std::map<std::string, modelTemplate1dType*>::iterator processEntry1d = processEntries1d_.begin();
	    processEntry1d != processEntries1d_.end(); ++processEntry1d ) {
	histograms.push_back(processEntry1d->second->histogram_);
	const RooAbsBinning& xRange = processEntry1d->second->xRef_->getBinning();
	double xMin = xRange.lowBound();
	double xMax = xRange.highBound();
	xRanges.push_back(double_pair(xMin, xMax));
      }

      std::string auxHistogramName = std::string(processName_).append("_").append("auxHistogram_1dPlus");
      delete auxHistogram_;
      auxHistogram_ = makeConcatenatedHistogram(auxHistogramName, histograms, xRanges);

//--- adjust range to be fitted for first variable
//    to x-axis of concatenated histogram
      const std::string varName = varNames_.front();
      RooRealVar* xRef_aux = processEntries1d_[varName]->xRef_;
      xRef_aux->setMin(auxHistogram_->GetXaxis()->GetXmin());
      xRef_aux->setMax(auxHistogram_->GetXaxis()->GetXmax());

      std::string auxDataHistName = std::string(processName_).append("_").append(varName).append("_rooDataHist");
      delete auxDataHist_;
      auxDataHist_ = new RooDataHist(auxDataHistName.data(), auxDataHistName.data(), *xRef_aux, auxHistogram_);

      std::string pdfName = std::string(processName_).append("_").append(varName).append("_rooHistPdf");
      pdf_ = new RooHistPdf(pdfName.data(), pdfName.data(), *xRef_aux, *auxDataHist_);
      
      pdfIsOwned_ = true;
    } else if ( fitMode_ == kNd ) {
      pdfName_ = std::string(processName_).append("_").append("pdf");

      TObjArray pdf1dCollection;
      for ( std::map<std::string, modelTemplate1dType*>::iterator processEntry1d = processEntries1d_.begin();
	    processEntry1d != processEntries1d_.end(); ++processEntry1d ) {
	pdf1dCollection.Add(processEntry1d->second->pdf1d_);
      }

      std::string pdfArgName = std::string(processName_).append("_pdfArgs");
      RooArgList pdfArgs(pdf1dCollection, pdfArgName.data());
      
      //std::cout << "--> creating RooProdPdf with name = " << pdfName_ << std::endl;
      pdf_ = new RooProdPdf(pdfName_.data(), pdfName_.data(), pdfArgs);

      pdfIsOwned_ = true;
    }
  }
}

void TemplateBgEstFit::modelTemplateNdType::fluctuate(bool fluctStat, bool fluctSys)
{
  for ( std::map<std::string, modelTemplate1dType*>::iterator processEntry1d = processEntries1d_.begin();
	processEntry1d != processEntries1d_.end(); ++processEntry1d ) {
    processEntry1d->second->fluctuate(fluctStat, fluctSys);
  }

  if ( pdfIsOwned_ ) delete pdf_;
  buildPdf();
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

TemplateBgEstFit::TemplateBgEstFit(const edm::ParameterSet& cfg)
  : error_(0)
{
  std::cout << "<TemplateBgEstFit::TemplateBgEstFit>:" << std::endl;

  edm::ParameterSet cfgFit = cfg.getParameter<edm::ParameterSet>("fit");

  std::string fitMode_string = ( cfgFit.exists("mode") ) ? cfgFit.getParameter<std::string>("mode") : fitMode_Nd;
  if ( fitMode_string == fitMode_1dPlus ) {
    fitMode_ = k1dPlus;
  } else if ( fitMode_string == fitMode_Nd ) {
    fitMode_ = kNd;
  } else {
    edm::LogError ("TemplateBgEstFit") << " Invalid 'fitMode' parameter = " << fitMode_string << " !!";
    error_ = 1;
  }

  //std::cout << " fitMode = " << fitMode_string << std::endl;

  std::map<std::string, bool> applyNormConstraints;
  std::map<std::string, double> meanNormConstraints;
  std::map<std::string, double> sigmaNormConstraints;
  if ( cfgFit.exists("constraints") ) {
    edm::ParameterSet cfgProcesses = cfgFit.getParameter<edm::ParameterSet>("constraints");
    vstring processNames = cfgProcesses.getParameterNamesForType<edm::ParameterSet>();
    for ( vstring::const_iterator processName = processNames.begin(); 
	  processName != processNames.end(); ++processName ) {
      edm::ParameterSet cfgProcess = cfgProcesses.getParameter<edm::ParameterSet>(*processName);
      
      applyNormConstraints[*processName] = true;
      meanNormConstraints[*processName] = cfgProcess.getParameter<double>("norm");
      sigmaNormConstraints[*processName] = cfgProcess.getParameter<double>("uncertainty");
    }
  }

  cutUnfittedRegion_ = ( cfgFit.exists("cutUnfittedRegion") ) ? cfgFit.getParameter<bool>("cutUnfittedRegion") : false;

  printLevel_ = ( cfgFit.exists("printLevel") ) ? cfgFit.getParameter<int>("printLevel") : 1;
  printWarnings_  = ( cfgFit.exists("printWarnings") ) ? cfgFit.getParameter<bool>("printWarnings") : true;

//--- read list of variables to be used in fit
//    (for each variable: name, title and range to be fitted)
  edm::ParameterSet cfgVariables = cfgFit.getParameter<edm::ParameterSet>("variables");
  vstring varNames = cfgVariables.getParameterNamesForType<edm::ParameterSet>();
  for ( vstring::const_iterator varName = varNames.begin(); 
	varName != varNames.end(); ++varName ) {
    edm::ParameterSet cfgVariable = cfgVariables.getParameter<edm::ParameterSet>(*varName);

    std::cout << "--> including variable = " << (*varName) << " in fit." << std::endl;
    varNames_.push_back(*varName);

    std::string name = cfgVariable.getParameter<std::string>("name");
    std::string title = cfgVariable.getParameter<std::string>("title");
    double xMin = cfgVariable.getParameter<double>("xMin");
    double xMax = cfgVariable.getParameter<double>("xMax");

    x_[*varName] = new RooRealVar(name.data(), title.data(), xMin, xMax);
  }

//--- read configuration parameters specifying signal and background templates
//    fitted to the distributions observed in (pseudo)data
//    (for each process: name; for each combination of variable and process: name of DQM MonitorElement holding template histogram)
  edm::ParameterSet cfgProcesses = cfg.getParameter<edm::ParameterSet>("processes");
  vstring processNames = cfgProcesses.getParameterNamesForType<edm::ParameterSet>();
  for ( vstring::const_iterator processName = processNames.begin(); 
	processName != processNames.end(); ++processName ) {
    edm::ParameterSet cfgProcess = cfgProcesses.getParameter<edm::ParameterSet>(*processName);

    processNames_.push_back(*processName);

    modelTemplateNdType* processEntry = new modelTemplateNdType(*processName, cfgProcess, fitMode_, cutUnfittedRegion_,
								applyNormConstraints[*processName], 
								meanNormConstraints[*processName], sigmaNormConstraints[*processName]);
    
    edm::ParameterSet cfgTemplates = cfgProcess.getParameter<edm::ParameterSet>("templates");
    for ( vstring::const_iterator varName = varNames_.begin();
	  varName != varNames_.end(); ++varName ) {
      if ( !cfgTemplates.exists(*varName) ) {
	edm::LogError ("TemplateBgEstFit") << " No Template of variable = " << (*varName) 
					   << " defined for process = " << (*processName) << " !!";
	error_ = 1;
	continue;
      }

      edm::ParameterSet cfgTemplate = cfgTemplates.getParameter<edm::ParameterSet>(*varName);

      std::string meName = cfgTemplate.getParameter<std::string>("meName");

      bool applySmoothing = cfgTemplate.exists("smoothing");
      edm::ParameterSet cfgSmoothing = ( applySmoothing ) ? 
	cfgTemplate.getParameter<edm::ParameterSet>("smoothing") : edm::ParameterSet();

      processEntry->addElement(*varName, x_[*varName], meName, applySmoothing, cfgSmoothing);
    }

    processEntries_[*processName] = processEntry;
  }
  
//--- read configuration parameters specifying distributions observed in (pseudo)data 
//    that are to be fitted
//    (for each variable: name of DQM MonitorElements holding template histogram)
//
//    WARNING: dataDistrNdType::addElement and modelTemplateNdType::addElement need to be called
//              with variable names sorted in the **exact** same order !!
//
  edm::ParameterSet cfgData = cfg.getParameter<edm::ParameterSet>("data");

  dataEntry_ = new dataDistrNdType(fitMode_, cutUnfittedRegion_);

  edm::ParameterSet cfgDistributions = cfgData.getParameter<edm::ParameterSet>("distributions");
  for ( vstring::const_iterator varName = varNames_.begin();
        varName != varNames_.end(); ++varName ) {
    if ( !cfgDistributions.exists(*varName) ) {
      edm::LogError ("TemplateBgEstFit") << " No Template of variable = " << (*varName) << " defined for data !!";
      error_ = 1;
      continue;
    }

    edm::ParameterSet cfgDistribution = cfgDistributions.getParameter<edm::ParameterSet>(*varName);

    std::string meName = cfgDistribution.getParameter<std::string>("meName");

    dataEntry_->addElement(*varName, x_[*varName], meName);
  }

//--- read configuration parameters specifying options for making control plots
  edm::ParameterSet cfgControlPlots = cfg.getParameter<edm::ParameterSet>("output").getParameter<edm::ParameterSet>("controlPlots");
  controlPlotsFileName_ = cfgControlPlots.getParameter<std::string>("fileName");

//--- read configuration parameters specifying how statistical and systematic uncertainties 
//    on normalization factors determined by fit get estimated
  edm::ParameterSet cfgStatErr = cfg.getParameter<edm::ParameterSet>("estStatUncertainties");
  statErrNumSamplings_ = cfgStatErr.getParameter<edm::ParameterSet>("numSamplings").getParameter<int>("stat");
  statErrChi2redMax_ = cfgStatErr.getParameter<double>("chi2redMax");  
  statErrPrintLevel_ = cfgStatErr.getParameter<edm::ParameterSet>("verbosity").getParameter<int>("printLevel");
  statErrPrintWarnings_ = cfgStatErr.getParameter<edm::ParameterSet>("verbosity").getParameter<bool>("printWarnings");

  edm::ParameterSet cfgSysErr = cfg.getParameter<edm::ParameterSet>("estSysUncertainties");
  edm::ParameterSet cfgSysFluctuations = cfgSysErr.getParameter<edm::ParameterSet>("fluctuations");
  vstring sysFluctNames = cfgSysFluctuations.getParameterNamesForType<edm::ParameterSet>();
  for ( vstring::const_iterator sysFluctName = sysFluctNames.begin(); 
	sysFluctName != sysFluctNames.end(); ++sysFluctName ) {
    edm::ParameterSet cfgSysFluct = cfgSysFluctuations.getParameter<edm::ParameterSet>(*sysFluctName);

    double pullRMS = cfgSysFluct.getParameter<double>("pullRMS");
    double pullMin = cfgSysFluct.getParameter<double>("pullMin");
    double pullMax = cfgSysFluct.getParameter<double>("pullMax");

    std::string fluctMode_string = cfgSysFluct.getParameter<std::string>("mode");
    int fluctMode_int = -1;
    if ( fluctMode_string == fluctMode_coherent ) {
      fluctMode_int = kCoherent;
    } else if ( fluctMode_string == fluctMode_incoherent ) {
      fluctMode_int = kIncoherent;
    } else {
      edm::LogError ("TemplateBgEstFit::TemplateBgEstFit") << " Invalid 'mode' parameter = " << fluctMode_string << " !!";
      error_ = 1;
    }

    edm::ParameterSet cfgProcesses = cfgSysFluct.getParameter<edm::ParameterSet>("meNames");

    for ( vstring::iterator processName = processNames_.begin();
	  processName != processNames_.end(); ++processName ) {
      if ( !cfgProcesses.exists(*processName) ) {
	edm::LogError ("TemplateBgEstFit") << " No Estimate of systematic Uncertainty = " << (*sysFluctName) 
					   << " defined for process = " << (*processName) << " !!";
	error_ = 1;
	continue;
      }

      edm::ParameterSet cfgProcess = cfgProcesses.getParameter<edm::ParameterSet>(*processName);

      for ( vstring::const_iterator varName = varNames_.begin();
	    varName != varNames_.end(); ++varName ) {
	if ( !cfgProcess.exists(*varName) ) {
	  edm::LogError ("TemplateBgEstFit") << " No Estimate of systematic Uncertainty = " << (*sysFluctName) 
					     << " on variable = " << (*varName) << " defined for process = " << (*processName) << " !!";
	  error_ = 1;
	}

	std::string meName = cfgProcess.getParameter<std::string>(*varName);

	sysFluctDefType sysFluctDef;
	sysFluctDef.fluctName_ = (*sysFluctName);
	sysFluctDef.meName_ = meName;
	sysFluctDef.pullRMS_ = pullRMS;
	sysFluctDef.pullMin_ = pullMin;
	sysFluctDef.pullMax_ = pullMax;
	sysFluctDef.fluctMode_ = fluctMode_int;

	processEntries_[*processName]->processEntries1d_[*varName]->sysErrFluctuations_.push_back(sysFluctDef);
      }
    }
  }
  sysErrNumStatSamplings_ = cfgSysErr.getParameter<edm::ParameterSet>("numSamplings").getParameter<int>("stat");
  sysErrNumSysSamplings_ = cfgSysErr.getParameter<edm::ParameterSet>("numSamplings").getParameter<int>("sys");
  sysErrChi2redMax_ = cfgSysErr.getParameter<double>("chi2redMax");  
  sysErrPrintLevel_ = cfgSysErr.getParameter<edm::ParameterSet>("verbosity").getParameter<int>("printLevel");
  sysErrPrintWarnings_ = cfgSysErr.getParameter<edm::ParameterSet>("verbosity").getParameter<bool>("printWarnings");

  fitModel_ = 0;
  fitResult_ = 0;
}

TemplateBgEstFit::~TemplateBgEstFit()
{
  std::cout << "<TemplateBgEstFit::~TemplateBgEstFit>:" << std::endl;
  for ( processEntryMap::iterator it = processEntries_.begin();
	it != processEntries_.end(); ++it ) {
    delete it->second;
  }

  delete dataEntry_;

  for ( realVarMap::iterator it = x_.begin();
	it != x_.end(); ++it ) {
    delete it->second;
  }

  delete fitModel_;
  delete fitResult_;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

void TemplateBgEstFit::buildFitModel()
{
  std::string fitModelName = "fitModel";
  
  TObjArray fitModel_pdfCollection;
  TObjArray fitModel_normCollection;
  for ( processEntryMap::iterator processEntry = processEntries_.begin();
	processEntry != processEntries_.end(); ++processEntry ) {
    fitModel_pdfCollection.Add(processEntry->second->pdf_);
    fitModel_normCollection.Add(processEntry->second->norm_);
  }

  std::string fitModel_pdfArgName = std::string("fitModel").append("_pdfArgs");
  RooArgList fitModel_pdfArgs(fitModel_pdfCollection, fitModel_pdfArgName.data());
  std::string fitModel_normArgName = std::string("fitModel").append("_normArgs");
  RooArgList fitModel_normArgs(fitModel_normCollection, fitModel_normArgName.data());

  //std::cout << "--> creating RooAddPdf with name = " << fitModelName << std::endl;
  fitModel_ = new RooAddPdf(fitModelName.data(), fitModelName.data(), fitModel_pdfArgs, fitModel_normArgs);
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

void TemplateBgEstFit::fit(bool saveFitResult, int printLevel, int printWarnings) 
{
//--- fit model to data

//--- build list of fit options
  RooLinkedList fitOptions;

//--- check if "external" constraints exist on normalization factors to be determined by fit
//    (specified by Gaussian probability density functions with mean and sigma obtained
//     e.g. by level of agreement between Monte Carlo simulation and number of events observed in background enriched samples)
  TObjArray normConstraints_pdfCollection;
  for ( processEntryMap::iterator processEntry = processEntries_.begin();
	processEntry != processEntries_.end(); ++processEntry ) {
    if ( processEntry->second->applyNormConstraint_ ) normConstraints_pdfCollection.Add(processEntry->second->pdfNormConstraint_);
  }

  if ( normConstraints_pdfCollection.GetEntries() > 0 ) {
    std::string normConstraints_pdfArgName = std::string("normConstraints").append("_pdfArgs");
    RooArgSet normConstraints_pdfArgs(normConstraints_pdfCollection, normConstraints_pdfArgName.data());
    
    // CMSSW_3_1_x only (because needs ROOT version 5.22)
    fitOptions.Add(new RooCmdArg(RooFit::ExternalConstraints(normConstraints_pdfArgs)));
  }

//--- check if results of fit are to be saved for later analysis
//    (only necessary when fitting for "central values", **not** when estimating statistical or systematic uncertainties)
  if ( saveFitResult ) fitOptions.Add(new RooCmdArg(RooFit::Save(true)));

//--- stop Minuit from printing lots of information 
//    about progress on fit and warnings
  fitOptions.Add(new RooCmdArg(RooFit::PrintLevel(printLevel)));
  fitOptions.Add(new RooCmdArg(RooFit::PrintEvalErrors(printWarnings)));
  fitOptions.Add(new RooCmdArg(RooFit::Warnings(printWarnings_)));

  RooFitResult* fitResult = fitModel_->fitTo(*dataEntry_->fitData_, fitOptions);
  if ( saveFitResult ) fitResult_ = fitResult;

//--- delete fit option objects
  int numCmdArgs = fitOptions.GetSize();
  for ( int iCmdArg = 0; iCmdArg < numCmdArgs; ++iCmdArg ) {
    delete fitOptions.At(iCmdArg);
  }
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

void TemplateBgEstFit::endJob()
{
  std::cout << "<TemplateBgEstFit::endJob>:" << std::endl;

//--- check that DQMStore service is available
  if ( !edm::Service<DQMStore>().isAvailable() ) {
    edm::LogError ("TemplateBgEstFit::endJob") << " Failed to access dqmStore" 
					       << " --> histogram will NOT be fitted !!";
    return;
  }

  //DQMStore& dqmStore = (*edm::Service<DQMStore>());
  //dqmStore.showDirStructure();

//--- check that configuration parameters contain no errors,
//    retrieve MonitorElements from DQMStore
//    and check that all DQM MonitorElements have successfully been retrieved
  for ( processEntryMap::iterator processEntry = processEntries_.begin();
	processEntry != processEntries_.end(); ++processEntry ) {
    processEntry->second->initialize();
    if ( processEntry->second->error_ ) error_ = 1;
  }

  dataEntry_->initialize();
  if ( dataEntry_->error_ ) error_ = 1;
  
  if ( error_ ) {
    edm::LogError ("TemplateBgEstFit::endJob") << " Error in Configuration ParameterSet" 
					       << " --> histogram will NOT be fitted !!";
    return;
  }
  
//--- configure RooFit structure;
//    print-out structure once configuration finished
  buildFitModel();

  std::cout << ">>> RootFit model used for generalized Matrix method Fit <<<" << std::endl;
  fitModel_->printCompactTree();

  std::cout << ">>> RootFit Parameters <<<" << std::endl;
  for ( vstring::const_iterator varName = varNames_.begin();
	varName != varNames_.end(); ++varName ) {
    std::cout << "for Variable = " << (*varName) << ":" << std::endl;
    fitModel_->getParameters(dataEntry_->dataEntries1d_[*varName]->dataHist_)->Print("v");
  }

  std::cout << ">>> RootFit Observables <<<" << std::endl;
  for ( vstring::const_iterator varName = varNames_.begin();
	varName != varNames_.end(); ++varName ) {
    std::cout << "for Variable = " << (*varName) << ":" << std::endl;
    fitModel_->getObservables(dataEntry_->dataEntries1d_[*varName]->dataHist_)->Print("v");
  }

//--- fit template shapes of different signal and background processes
//    to distribution observed in (pseudo)data
  fit(true, printLevel_, printWarnings_);
  assert(fitResult_ != NULL);

//--- print-out fit results
  std::cout << ">>> Fit Results <<<" << std::endl;
  std::cout << " fitStatus = " << fitResult_->status() << std::endl;
  std::cout << " Chi2red = " << compChi2red() << std::endl;
  print(std::cout);

//--- save normalization factors determined and covariance matrix estimated by fit
  unpackFitResult(fitResult_, fitResultMean_, fitResultCov_);

//--- produce plot of different signal and background processes
//    using scale factors determined by fit
//    compared to distribution of (pseudo)data
  if ( controlPlotsFileName_ != "" ) makeControlPlots();

//--- estimate statistical uncertainties
  std::cout << ">>> Statistical Uncertainties <<<" << std::endl;
  estimateUncertainties(true, statErrNumSamplings_, false, 1, statErrChi2redMax_, 
			"estStatUncertainties", statErrPrintLevel_, statErrPrintWarnings_);
 
//--- estimate systematic uncertainties
  std::cout << ">>> Systematic Uncertainties <<<" << std::endl;
  estimateUncertainties(false, sysErrNumStatSamplings_, true, sysErrNumSysSamplings_, sysErrChi2redMax_,
			"estSysUncertainties", sysErrPrintLevel_, sysErrPrintWarnings_);

//--- estimate total (statistical + systematic) uncertainties
  std::cout << ">>> Total (statistical + systematic) Uncertainties <<<" << std::endl;
  double chi2redMax = TMath::Max(statErrChi2redMax_, sysErrChi2redMax_);
  int totErrPrintLevel = TMath::Min(statErrPrintLevel_, sysErrPrintLevel_);
  bool totErrPrintWarnings = (statErrPrintWarnings_ && sysErrPrintWarnings_);
  estimateUncertainties(true, sysErrNumStatSamplings_, true, sysErrNumSysSamplings_, chi2redMax,
			"estTotUncertainties", totErrPrintLevel, totErrPrintWarnings);
  
  std::cout << "done." << std::endl;
}

void TemplateBgEstFit::print(std::ostream& stream)
{
  stream << "Fit Parameter:" << std::endl;
  for ( processEntryMap::iterator processEntry = processEntries_.begin();
	processEntry != processEntries_.end(); ++processEntry ) {
    const std::string& processName = processEntry->first;
    RooRealVar* processNorm = processEntry->second->norm_;

//--- correct normalization factor determined by fit
//    for events events passing final analysis selection criteria
//    that are outside fitted region
    double normCorrFactor = dataEntry_->normCorrFactor_;

    stream << " " << processName << ": normalization = " << normCorrFactor*processNorm->getVal()
	   << " (within fitted region = " << processNorm->getVal() << ")";

    if ( processNorm->hasAsymError() ) {
      stream << " + " << normCorrFactor*processNorm->getAsymErrorHi() << " (" << processNorm->getAsymErrorHi() << ")"
	     << " - " << fabs(normCorrFactor*processNorm->getAsymErrorLo()) << " (" << fabs(processNorm->getAsymErrorLo()) << ")";
    } else if ( processNorm->hasError() ) {
      stream << " +/- " << normCorrFactor*processNorm->getError() << " (" << processNorm->getError() << ")";
    }

    stream << std::endl;
  }
}

void TemplateBgEstFit::makeControlPlots()
{
//--- produce control plot of distribution observed in (pseudo)data
//    versus sum of signal and background templates using normalization determined by fit

//--- stop ROOT from opening X-window for canvas output
//    (in order to be able to run in batch mode) 
  gROOT->SetBatch(true);

//--- produce control plots of smoothing functions
//    fitted to shape templates in order to reduce effect of statistical fluctuations
//    of shape templates on fit results
//    (and in particular effect of different statistical precision with which shape templates
//     are determined for different background processes in background enriched samples)
  //makeControlPlotsSmoothing();

//--- produce control plots of one and two sigma error contours 
//    showing correlation of estimated normalization factors
  const RooArgList& fitParameter = fitResult_->floatParsFinal();
  int numFitParameter = fitParameter.getSize();
  vstring labels(numFitParameter);
  for ( int iX = 0; iX < numFitParameter; ++iX ) {
    const RooAbsArg* paramX_arg = fitParameter.at(iX);
    labels[iX] = paramX_arg->GetName();
  }

  makeControlPlotsCovariance(fitResultMean_, fitResultMean_, fitResultCov_, labels, controlPlotsFileName_, "");

//--- produce control plots of distributions observed in (pseudo)data
//    compared to sum of signal and background templates
//    scaled by normalization factors determined by the fit
  makeControlPlotsObsDistribution();
}

void TemplateBgEstFit::makeControlPlotsSmoothing()
{
  TCanvas canvas("TemplateBgEstFit", "TemplateBgEstFit", defaultCanvasSizeX, defaultCanvasSizeY);
  canvas.SetFillColor(10);

  int defStyle_optStat = gStyle->GetOptStat();
  int defStyle_optFit = gStyle->GetOptStat();
  float defStyle_labelSizeX = gStyle->GetLabelSize("x");
  float defStyle_labelSizeY = gStyle->GetLabelSize("y");

  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(111);
  gStyle->SetLabelSize(0.03, "x");
  gStyle->SetLabelSize(0.03, "y");
  
  for ( processEntryMap::iterator processEntry = processEntries_.begin();
	processEntry != processEntries_.end(); ++processEntry ) {
    const std::string& processName = processEntry->first;
    
    for ( vstring::const_iterator varName = varNames_.begin(); 
	  varName != varNames_.end(); ++varName ) {

      modelTemplate1dType* processEntry1d = processEntry->second->processEntries1d_[*varName];

      if ( !processEntry1d->applySmoothing_ ) continue;

      std::string histogramName = std::string(processEntry1d->fluctHistogram_->GetName()).append("_cloned");
      TH1* histogram_cloned = (TH1*)processEntry1d->fluctHistogram_->Clone(histogramName.data());
      histogram_cloned->GetXaxis()->SetTitle(varName->data());
      histogram_cloned->SetLineStyle(1);
      histogram_cloned->SetLineColor(1);
      histogram_cloned->SetLineWidth(2);

      std::string tf1Name = std::string(processEntry1d->auxTF1Wrapper_->getTF1()->GetName()).append("_cloned");
      TF1* tf1_cloned = (TF1*)processEntry1d->auxTF1Wrapper_->getTF1()->Clone(tf1Name.data());
      tf1_cloned->SetLineStyle(1);
      tf1_cloned->SetLineColor(2);
      tf1_cloned->SetLineWidth(2);
      
      double yMax = histogram_cloned->GetMaximum();
      histogram_cloned->SetMaximum(1.3*yMax);

      TLegend legend(0.63, 0.38, 0.89, 0.54);
      legend.SetBorderSize(0);
      legend.SetFillColor(0);

      legend.AddEntry(histogram_cloned, "Template Histogram", "p");
      legend.AddEntry(tf1_cloned, "Smoothing Function", "l");

      histogram_cloned->Draw("e1p");
      tf1_cloned->Draw("lsame");
      
      //legend.Draw();

      canvas.Update();

      int errorFlag = 0;
      std::string fileNameParam = std::string(processName).append("_").append(*varName).append("_smoothing");
      std::string fileName = replace_string(controlPlotsFileName_, plotKeyword, fileNameParam, 1, 1, errorFlag);
      if ( !errorFlag ) {
	canvas.Print(fileName.data());
      } else {
	edm::LogError("TemplateBgEstFit::makeControlPlotsSmoothing") << " Failed to decode controlPlotsFileName = " 
								     << controlPlotsFileName_ << " --> skipping !!";
	return;
      }

      delete histogram_cloned;
      delete tf1_cloned;
    }
  }

  gStyle->SetOptStat(defStyle_optStat);
  gStyle->SetOptFit(defStyle_optFit);
  gStyle->SetLabelSize(defStyle_labelSizeX, "x");
  gStyle->SetLabelSize(defStyle_labelSizeY, "y");
}

void drawErrorEllipse(double x0, double y0, double errX0, double errY0, double Sxx, double Sxy, double Syy, 
		      const char* labelX, const char* labelY, const char* fileName)
{
//--- draw best fit estimate (x0,y0) 
//    plus one and two sigma error contours centered at (errX0,errY0) 
//    and with (correlated) uncertainties estimated by elements Sxx, Sxy, Syy of covariance matrix 
//    passed as function arguments
//    (note that since the covariance matrix is symmetric, 
//     there is no need to pass element Syx of the covariance matrix)

  TCanvas canvas("drawErrorEllipse", "drawErrorEllipse", 600, 600);
  canvas.SetFillColor(10);

//--- compute angle between first principal axis of error ellipse
//    and x-axis
  double alpha = 0.5*TMath::ATan2(-2*Sxy, Sxx - Syy);

  //std::cout << "alpha = " << alpha*180./TMath::Pi() << std::endl;

  double sinAlpha = TMath::Sin(alpha);
  double cosAlpha = TMath::Cos(alpha);

//--- compute covariance axis in coordinate system
//    defined by principal axes of error ellipse
  double Suu = Sxx*sinAlpha*sinAlpha + 2*Sxy*sinAlpha*cosAlpha + Syy*cosAlpha*cosAlpha;
  double Svv = Sxx*cosAlpha*cosAlpha + 2*Sxy*sinAlpha*cosAlpha + Syy*sinAlpha*sinAlpha;
  
//--- resolve ambiguity which axis represents the first principal axis
//    and which represents the second principal axis
//
//    NOTE: in case Sxy > 0. (correlation of variables X and Y), 
//          the principal axis needs to point in direction of either the first or the third quadrant;
//          in case Sxy < 0. (anti-correlation of variables X and Y), 
//          the principal axis needs to point in direction of either the second or the fourth quadrant.
  double sigmaX_transformed = 0.;
  double sigmaY_transformed = 0.;
  if ( (Sxy >= 0. && TMath::Abs(alpha) <= 0.5*TMath::Pi()) || 
       (Sxy <  0. && TMath::Abs(alpha) >  0.5*TMath::Pi()) ) {
    sigmaX_transformed = TMath::Sqrt(TMath::Max(Suu, Svv));
    sigmaY_transformed = TMath::Sqrt(TMath::Min(Suu, Svv));
  } else {
    sigmaX_transformed = TMath::Sqrt(TMath::Min(Suu, Svv));
    sigmaY_transformed = TMath::Sqrt(TMath::Max(Suu, Svv));
  }

  TEllipse oneSigmaErrorEllipse(errX0, errY0, sigmaX_transformed*1., sigmaY_transformed*1., 0., 360., alpha*180./TMath::Pi()); 
  oneSigmaErrorEllipse.SetFillColor(5);
  oneSigmaErrorEllipse.SetLineColor(44);
  oneSigmaErrorEllipse.SetLineWidth(1);
  TEllipse twoSigmaErrorEllipse(errX0, errY0, sigmaX_transformed*2., sigmaY_transformed*2., 0., 360., alpha*180./TMath::Pi()); 
  TSeqCollection* colors = gROOT->GetListOfColors();
  if ( colors && colors->At(42) ) {
    TColor* orange = (TColor*)colors->At(42);
    orange->SetRGB(1.00,0.80,0.00);
  } else {
    edm::LogWarning ("drawErrorEllipse") << " Failed to access list of Colors from gROOT object"
					 << " --> skipping definition of Color 'orange' !!";
  }
  twoSigmaErrorEllipse.SetFillColor(42);
  twoSigmaErrorEllipse.SetLineColor(44);
  twoSigmaErrorEllipse.SetLineWidth(1);

  TMarker centralValueMarker(x0, y0, 5);
  centralValueMarker.SetMarkerSize(2);

//--- create dummy histogram  
//    defining region to be plotted
  double minX = TMath::Min(errX0 - 2.2*(TMath::Abs(cosAlpha*sigmaX_transformed) + TMath::Abs(sinAlpha*sigmaY_transformed)),
			   x0 - 0.2*(TMath::Abs(cosAlpha*sigmaX_transformed) + TMath::Abs(sinAlpha*sigmaY_transformed)));
  double maxX = TMath::Max(errX0 + 2.8*(TMath::Abs(cosAlpha*sigmaX_transformed) + TMath::Abs(sinAlpha*sigmaY_transformed)),
			   x0 + 0.2*(TMath::Abs(cosAlpha*sigmaX_transformed) + TMath::Abs(sinAlpha*sigmaY_transformed)));
  double minY = TMath::Min(errY0 - 2.2*(TMath::Abs(sinAlpha*sigmaX_transformed) + TMath::Abs(cosAlpha*sigmaY_transformed)),
			   y0 - 0.2*(TMath::Abs(sinAlpha*sigmaX_transformed) + TMath::Abs(cosAlpha*sigmaY_transformed)));
  double maxY = TMath::Max(errY0 + 2.8*(TMath::Abs(sinAlpha*sigmaX_transformed) + TMath::Abs(cosAlpha*sigmaY_transformed)),
			   y0 + 0.2*(TMath::Abs(sinAlpha*sigmaX_transformed) + TMath::Abs(cosAlpha*sigmaY_transformed)));
			   
  if ( TMath::Abs(maxX - minX) < epsilon || 
       TMath::Abs(maxY - minY) < epsilon ) {
    if ( TMath::Abs(maxX - minX) < epsilon ) edm::LogWarning ("drawErrorEllipse") << " Invalid x-range: minX = maxX = " << minX;
    if ( TMath::Abs(maxY - minY) < epsilon ) edm::LogWarning ("drawErrorEllipse") << " Invalid y-range: minY = maxY = " << minY;
    edm::LogWarning ("drawErrorEllipse") << " --> skipping drawing of Error ellipse for labelX = " << labelX << ","
					 << " labelY = " << labelY << " !!";
    return;
  }

//--- create dummy histogram  
  TH2F dummyHistogram("dummyHistogram", "dummyHistogram", 5, minX, maxX, 5, minY, maxY);
  dummyHistogram.SetTitle("");
  dummyHistogram.SetStats(false);
  dummyHistogram.SetXTitle(labelX);
  dummyHistogram.SetYTitle(labelY);
  dummyHistogram.SetTitleOffset(1.35, "Y");

  dummyHistogram.Draw("AXIS");
  
  twoSigmaErrorEllipse.Draw();
  oneSigmaErrorEllipse.Draw();

  centralValueMarker.Draw();

  TLegend legend(0.70, 0.70, 0.89, 0.89);
  legend.SetBorderSize(0);
  legend.SetFillColor(0);
  legend.AddEntry(&centralValueMarker, "best Fit value", "p");
  legend.AddEntry(&oneSigmaErrorEllipse, "1#sigma Contour", "f");
  legend.AddEntry(&twoSigmaErrorEllipse, "2#sigma Contour", "f");

  legend.Draw();

  canvas.Print(fileName);
}

void TemplateBgEstFit::makeControlPlotsCovariance(TVectorD estimate, TVectorD errMean, TMatrixD errCov, const vstring& labels, 
						  const std::string& controlPlotsFileName, const char* type)
{
  int numFitParameter = estimate.GetNoElements();
  assert(errMean.GetNoElements() == numFitParameter);
  assert(errCov.GetNrows() == numFitParameter && errCov.GetNcols() == numFitParameter);
  for ( int iX = 0; iX < numFitParameter; ++iX ) {
    double x0 = estimate(iX);
    double errX0 = errMean(iX);
    double Sxx = errCov(iX, iX);
    const char* labelX = labels[iX].data();

    for ( int iY = 0; iY < iX; ++iY ) {
      double y0 = estimate(iY);
      double errY0 = errMean(iY);
      double Syy = errCov(iY, iY);
      const char* labelY = labels[iY].data();

      double Sxy = errCov(iX, iY);
      std::string fileNameParam = std::string("corr_").append(labelX).append("_vs_").append(labelY);
      if ( type != "" ) fileNameParam.append("_").append(type);
      
      int errorFlag = 0;
      std::string fileName = replace_string(controlPlotsFileName, plotKeyword, fileNameParam, 1, 1, errorFlag);
      if ( !errorFlag ) {
	drawErrorEllipse(x0, y0, errX0, errY0, Sxx, Sxy, Syy, labelX, labelY, fileName.data());
      } else {
	edm::LogError("drawErrorEllipses") << " Failed to decode controlPlotsFileName = " << controlPlotsFileName 
					   << " --> skipping !!";
	return;
      }
    }
  }
}

void TemplateBgEstFit::makeControlPlotsObsDistribution()
{
  TCanvas canvas("TemplateBgEstFit", "TemplateBgEstFit", defaultCanvasSizeX, defaultCanvasSizeY);
  canvas.SetFillColor(10);

  for ( vstring::const_iterator varName = varNames_.begin(); 
	varName != varNames_.end(); ++varName ) {
    TH1* fittedHistogram_sum = 0;

    std::vector<TH1*> fittedHistograms;

    TLegend legend(0.67, 0.63, 0.89, 0.89);
    legend.SetBorderSize(0);
    legend.SetFillColor(0);

    for ( processEntryMap::iterator processEntry = processEntries_.begin();
	  processEntry != processEntries_.end(); ++processEntry ) {
      const std::string& processName = processEntry->first;
      double processNorm = processEntry->second->norm_->getVal();

//--- correct normalization factor determined by fit
//    for events events passing final analysis selection criteria
//    that are outside fitted region
      double normCorrFactor = dataEntry_->normCorrFactor_;

      TH1* processShapeHist = processEntry->second->processEntries1d_[*varName]->me_->getTH1();

      std::string fittedHistogramName_process = std::string(processShapeHist->GetName()).append("_cloned");
      TH1* fittedHistogram_process = (TH1*)processShapeHist->Clone(fittedHistogramName_process.data());

      if ( getIntegral(fittedHistogram_process) > 0. ) {
	fittedHistogram_process->Scale(normCorrFactor*processNorm/getIntegral(fittedHistogram_process));
      }

      fittedHistogram_process->SetLineColor(processEntry->second->lineColor_);
      fittedHistogram_process->SetLineStyle(processEntry->second->lineStyle_);
      fittedHistogram_process->SetLineWidth(processEntry->second->lineWidth_);
      
      fittedHistograms.push_back(fittedHistogram_process);

      legend.AddEntry(fittedHistogram_process, processName.data(), "l");

      if ( !fittedHistogram_sum ) {
	std::string fittedHistogramName_sum = std::string(*varName).append("_histogram_fittedSum");
	fittedHistogram_sum = (TH1*)fittedHistogram_process->Clone(fittedHistogramName_sum.data());
	fittedHistogram_sum->SetStats(false);
	fittedHistogram_sum->GetXaxis()->SetTitle(varName->data());
	fittedHistogram_sum->SetLineColor(1); // black
	fittedHistogram_sum->SetLineStyle(1); // solid
	fittedHistogram_sum->SetLineWidth(fittedHistogram_process->GetLineWidth());
      } else {
	fittedHistogram_sum->Add(fittedHistogram_process);
      }
    }

    TH1* dataHistogram = dataEntry_->dataEntries1d_[*varName]->me_->getTH1();
    dataHistogram->SetMarkerStyle(8);
    legend.AddEntry(dataHistogram, "final Evt. Sel.", "p");

    double yMax = TMath::Max(fittedHistogram_sum->GetMaximum(),
			     dataHistogram->GetMaximum());
    fittedHistogram_sum->SetMaximum(1.3*yMax);
    legend.AddEntry(fittedHistogram_sum, "fitted #Sigma", "l");

    fittedHistogram_sum->Draw("hist");

    for ( std::vector<TH1*>::const_iterator fittedHistogram_process = fittedHistograms.begin();
	  fittedHistogram_process != fittedHistograms.end(); ++fittedHistogram_process ) {
      (*fittedHistogram_process)->Draw("histsame");
    }

    dataHistogram->Draw("e1psame");

    legend.Draw();
  
    canvas.Update();

    int errorFlag = 0;
    std::string fileName = replace_string(controlPlotsFileName_, plotKeyword, *varName, 1, 1, errorFlag);
    if ( !errorFlag ) {
      canvas.Print(fileName.data());
    } else {
      edm::LogError("TemplateBgEstFit::makeControlPlotsObsDistribution") << " Failed to decode controlPlotsFileName = " 
									 << controlPlotsFileName_ << " --> skipping !!";
      return;
    }

    delete fittedHistogram_sum;
    for ( std::vector<TH1*>::iterator it = fittedHistograms.begin();
	  it != fittedHistograms.end(); ++it ) {
      delete (*it);
    }
  }
}

double TemplateBgEstFit::compChi2red() 
{
  //std::cout << "<TemplateBgEstFit::compChi2red>:" << std::endl;
  
  double chi2 = 0.;
  int numDoF = 0;

  for ( vstring::const_iterator varName = varNames_.begin(); 
	varName != varNames_.end(); ++varName ) {
    const TH1* histogramData = dataEntry_->dataEntries1d_[*varName]->me_->getTH1();
    const RooAbsBinning& xRange = x_[*varName]->getBinning();
    double xMin = xRange.lowBound();
    double xMax = xRange.highBound();
    int numBins = histogramData->GetNbinsX();
    for ( int iBin = 0; iBin < (numBins + 2); ++iBin ) {
      double dataBinCenter = histogramData->GetBinCenter(iBin);

//--- restrict computation of chi^2 to region included in fit
      if ( !(dataBinCenter > xMin && dataBinCenter < xMax) ) continue;

      double dataBinContent = histogramData->GetBinContent(iBin);
      double dataBinError = histogramData->GetBinError(iBin);
      
      double fitBinContent = 0.;
      double fitBinError2 = 0.;
      for ( processEntryMap::iterator processEntry = processEntries_.begin();
	    processEntry != processEntries_.end(); ++processEntry ) {
	const TH1* histogramProcess = processEntry->second->processEntries1d_[*varName]->me_->getTH1();

	double processBinContent = histogramProcess->GetBinContent(iBin);
	double processBinError = histogramProcess->GetBinError(iBin);

	double processNorm = processEntry->second->norm_->getVal();
	double processNormCorrFactor = processEntry->second->normCorrFactor_;
	//std::cout << "processName = " << processEntry->first << ": processNormCorrFactor = " << processNormCorrFactor << std::endl;

	double processBinContent_scaled = processNormCorrFactor*processNorm*processBinContent;
	double processBinError_scaled = processNormCorrFactor*processNorm*processBinError;

	fitBinContent += processBinContent_scaled;
	fitBinError2 += processBinError_scaled*processBinError_scaled;
      }
      
      //std::cout << "iBin = " << iBin << ": dataBinContent = " << dataBinContent << ", fitBinContent = " << fitBinContent << std::endl;

      double diffBinContent2 = (dataBinContent - fitBinContent)*(dataBinContent - fitBinContent);
      double diffBinError2 = fitBinError2 + dataBinError*dataBinError;
      
      if ( diffBinError2 > 0. ) {
	chi2 += (diffBinContent2/diffBinError2);
	++numDoF;
      }
    }
  }

//--- correct number of degrees of freedom
//    for number of fitted parameters
  numDoF -= processEntries_.size();

  //std::cout << "chi2 = " << chi2 << std::endl;
  //std::cout << "numDoF = " << numDoF << std::endl;
  
  if ( numDoF > 0 ) {
    return (chi2/numDoF);
  } else {
    edm::LogWarning ("compChi2red") << " numDoF = " << numDoF << " must not be negative"
				    << " returning Chi2red = 1.e+3 !!";
    return 1.e+3;
  }
}

void TemplateBgEstFit::estimateUncertainties(bool fluctStat, int numStatSamplings, bool fluctSys, int numSysSamplings, 
					     double chi2redMax, const char* type, int printLevel, bool printWarnings)
{
  int numProcesses = processEntries_.size();
  TVectorD fitValues(numProcesses);
  BgEstMean mean(numProcesses);
  BgEstMedian median(numProcesses);
  BgEstCovMatrix cov(numProcesses);

  unsigned numTotFits = 0;
  unsigned numGoodFits = 0;

  for ( int iRndStat = 0; iRndStat < numStatSamplings; ++iRndStat ) {
    for ( int iRndSys = 0; iRndSys < numSysSamplings; ++iRndSys ) {

      std::cout << "<TemplateBgEstFit::estimateUncertainties>: iRndStat = " << iRndStat << ", iRndSys = " << iRndSys << std::endl;

//--- fluctuate distributions observed in (pseudo)data
      dataEntry_->fluctuate(true, false);

//--- fluctuate template histograms fitted to the (pseudo)data  
//
//    CV: treat statistical uncertainties on template histograms
//        as systematic uncertainties of background estimation
//
      for ( processEntryMap::iterator processEntry = processEntries_.begin();
            processEntry != processEntries_.end(); ++processEntry ) {
        //processEntry->second->fluctuate(fluctStat, fluctSys);
	processEntry->second->fluctuate(fluctSys, fluctSys);
      }

      delete fitModel_;
      buildFitModel();

      fit(false, printLevel, printWarnings);
      ++numTotFits;

      for ( int iProcess = 0; iProcess < numProcesses; ++iProcess ) {
	const std::string& processName = processNames_[iProcess];
	fitValues(iProcess) = processEntries_[processName]->norm_->getVal();
	//std::cout << " fitValue(iProcess = " << iProcess << ", processName = " << processName << ")"
	//	    << " = " << fitValues(iProcess) << std::endl;
      }

      double chi2red = compChi2red();
      //std::cout << "Chi2red = " << chi2red << std::endl;
      if ( !(chi2red < chi2redMax) ) continue;

      mean.update(fitValues);
      median.update(fitValues);
      cov.update(fitValues);

      ++numGoodFits;
    }
  }

  double badFitFraction = (numTotFits - numGoodFits)/((double)numTotFits);
  std::cout << "fraction of Samplings discarded due to bad Fit quality = " << badFitFraction << std::endl;

  std::cout << "Mean:" << std::endl;
  mean.print(std::cout, &processNames_);
  std::cout << "Median:" << std::endl;
  median.print(std::cout, &processNames_);
  std::cout << "Covariance Matrix:" << std::endl;
  cov.print(std::cout, &processNames_);
  
  if ( controlPlotsFileName_ != "" ) {
    //makeControlPlotsCovariance(fitResultMean_, median(), cov(), processNames_, controlPlotsFileName_, type);
    makeControlPlotsCovariance(fitResultMean_, fitResultMean_, cov(), processNames_, controlPlotsFileName_, type);
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(TemplateBgEstFit);
