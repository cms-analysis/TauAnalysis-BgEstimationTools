#include "TauAnalysis/BgEstimationTools/plugins/GenMatrixBgEstFit.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <TObjArray.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TROOT.h>

#include <RooRealVar.h>
#include <RooConstVar.h>
#include <RooFormulaVar.h>
#include <RooProdPdf.h>
#include <RooAddPdf.h>
#include <RooArgSet.h>
#include <RooGlobalFunc.h>
#include <RooPlot.h>
#include <RooFit.h>

#include <iomanip>
#include <fstream>

const int defaultCanvasSizeX = 800;
const int defaultCanvasSizeY = 600;

const float maxNorm = 1.e+6;

float getRandom()
{
  static TRandom3 u;

  float value = -10.;
  while ( value < -3 || value > +3 ) {
    value = u.Rndm();
  }

  return value;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

void getProbability(double probDensity1, double probDensity2, const Float_t* binBoundaries,
		    double& prob1, double& prob2)
{
//--- "probability" values used by GenMatrixBgEst1dPdf are in fact probability densities,
//    so need to integrate probability density over bin-width 

  //std::cout << "<getProbability>:" << std::endl;
  //std::cout << " probDensity1 = " << probDensity1 << std::endl;
  //std::cout << " probDensity2 = " << probDensity2 << std::endl;

  double prob1_unnormalized = (binBoundaries[1] - binBoundaries[0])*probDensity1;
  double prob2_unnormalized = (binBoundaries[2] - binBoundaries[1])*probDensity2;

  //std::cout << " prob1_unnormalized = " << prob1_unnormalized << std::endl;
  //std::cout << " prob2_unnormalized = " << prob2_unnormalized << std::endl;
  
  prob1 = prob1_unnormalized/(prob1_unnormalized + prob2_unnormalized);
  prob2 = prob2_unnormalized/(prob1_unnormalized + prob2_unnormalized);

  //std::cout << " prob1 = " << prob1 << std::endl;
  //std::cout << " prob2 = " << prob2 << std::endl;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

GenMatrixBgEstFit::objVarEntryType::objVarEntryType(const std::string& name, float xBoundary, float xMin, float xMax)
  : name_(name), xBoundary_(xBoundary), xMin_(xMin), xMax_(xMax) 
{ 
  //std::cout << "<objVarEntryType::objVarEntryType>:" << std::endl;
  //std::cout << " name = " << name_ << std::endl;
  //std::cout << " xMin = " << xMin_ << std::endl;
  //std::cout << " xMax = " << xMax_ << std::endl;

  std::string xName = std::string(name_);
  std::string xTitle = xName;
  std::cout << "--> creating RooRealVar with name = " << xName << ","
	    << " range of validity = [" << xMin_ << ", " << xMax_ << "]..." << std::endl;
  x_ = new RooRealVar(xName.data(), xTitle.data(), xMin, xMax);
}

GenMatrixBgEstFit::objVarEntryType::~objVarEntryType()
{
  //std::cout << "<objVarEntryType::~objVarEntryType>:" << std::endl;

  delete x_;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

GenMatrixBgEstFit::pdf1dEntryType::pdf1dEntryType(const std::string& name, 
						  objVarEntryType* objVarEntry, 
						  float prob, bool fixProb)
  : name_(name), fixProb_(fixProb)
{
  //std::cout << "<pdf1dEntryType::pdf1dEntryType>:" << std::endl;
  //std::cout << " name = " << name_ << std::endl;
  
  const Int_t numBins = 2;  
  Float_t binBoundaries[numBins + 1];
  binBoundaries[0] = objVarEntry->xMin_;
  binBoundaries[1] = objVarEntry->xBoundary_;
  binBoundaries[2] = objVarEntry->xMax_;

//--- "probability" values used by GenMatrixBgEst1dPdf are in fact probability densities,
//    so need to devide probability by bin-width (and renormalize)
//
//    NOTE: from the condition 
//
//         1 = prob + (1 - prob) = prob1 + prob2
//           = (binBoundaries[1] - binBoundaries[0])*probDensity1 + (binBoundaries[2] - binBoundaries[1])*probDensity2 
// 
//          it follows that
//
//         probDensity2 =                 1.                    binBoundaries[1] - binBoundaries[0] 
//                        ----------------------------------- - ----------------------------------- * probDensity1 
//                        binBoundaries[2] - binBoundaries[1]   binBoundaries[2] - binBoundaries[1]
//
  double probDensity1 = prob/(binBoundaries[1] - binBoundaries[0]);
  double probDensity2_c0 = 1./(binBoundaries[2] - binBoundaries[1]);
  double probDensity2_c1 = (binBoundaries[1] - binBoundaries[0])/(binBoundaries[2] - binBoundaries[1]);

  //std::cout << " probDensity1*(binBoundaries[1] - binBoundaries[0]) = " 
  //	      << probDensity1*(binBoundaries[1] - binBoundaries[0]) << std::endl;
  //double probDensity2 = probDensity2_c0 - probDensity2_c1*probDensity1;
  //std::cout << " probDensity2*(binBoundaries[2] - binBoundaries[1]) = " 
  //	      << probDensity2*(binBoundaries[2] - binBoundaries[1]) << std::endl;

  std::string prob1Name = std::string(name_).append("_prob1");
  std::string prob1Title = prob1Name;
  if ( fixProb ) {
    std::cout << "--> creating RooConstVar with name = " << prob1Name << ":" 
	      << " probability = " << prob << ", prob. density = " << probDensity1 << "..." << std::endl;
    prob1_ = new RooConstVar(prob1Name.data(), prob1Title.data(), probDensity1);
    //std::cout << "--> creating RooConstVar with name = " << prob1Name << ":" 
    //	        << " probability = " << prob << "..." << std::endl;
    //prob1_ = new RooConstVar(prob1Name.data(), prob1Title.data(), prob);
  } else {
    std::cout << "--> creating RooRealVar with name = " << prob1Name << ", start-values:" 
	      << " probability = " << prob << ", prob. density = " << probDensity1 << "..." << std::endl;
    prob1_ = new RooRealVar(prob1Name.data(), prob1Title.data(), probDensity1, 0., 1./(binBoundaries[1] - binBoundaries[0]));
    //std::cout << "--> creating RooRealVar with name = " << prob1Name << ":" 
    //	        << " probability start-value = " << prob << "..." << std::endl;
    //prob1_ = new RooRealVar(prob1Name.data(), prob1Title.data(), prob, 0., 1.);
  }

  std::string prob2Name = std::string(name_).append("_prob2");
  std::string prob2Title = prob2Name;
  std::ostringstream prob2Function;
  prob2Function << probDensity2_c0 << "-" << "(" << probDensity2_c1 << "*" << prob1Name << ")";
  //prob2Function << "1-" << prob1Name;
  //std::cout << " prob2Function = " << prob2Function.str() << std::endl;
  std::cout << "--> creating RooFormulaVar with name = " << prob2Name << "," 
	    << " depending on variable = " << prob1_->GetName() 
	    << " via function = " << prob2Function.str() << "..." << std::endl;
  prob2_ = new RooFormulaVar(prob2Name.data(), prob2Title.data(), prob2Function.str().data(), RooArgList(*prob1_));

  std::string pdf1dName = std::string(name_).append("_pdf1d");
  std::string pdf1dTitle = pdf1dName;
  TObjArray binProbabilities;
  binProbabilities.Add(prob1_);
  binProbabilities.Add(prob2_);
  std::cout << "--> creating GenMatrixBgEst1dPdf with name = " << pdf1dName.data() << ","
	    << " depending on variables = " << prob1_->GetName() << ", " << prob2_->GetName() << std::endl;
  pdf1d_ = new GenMatrixBgEst1dPdf(pdf1dName.data(), pdf1dTitle.data(), *objVarEntry->x_, numBins, binBoundaries, binProbabilities);
}

GenMatrixBgEstFit::pdf1dEntryType::~pdf1dEntryType()
{
  //std::cout << "<pdf1dEntryType::~pdf1dEntryType>:" << std::endl;

  delete pdf1d_;

  delete prob1_;
  delete prob2_;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

GenMatrixBgEstFit::pdfSingleProcessEntryType::pdfSingleProcessEntryType(const std::string& name, 
									const std::vector<objVarEntryType*>& objVarEntries,
									const std::vector<float>& prob, const std::vector<bool>& fixProb,
									float norm, bool fixNorm)
  : name_(name),
    lineColor_(0), lineStyle_(0), lineWidth_(0)
{
  //std::cout << "<pdfSingleProcessEntryType::pdfSingleProcessEntryType>:" << std::endl;
  //std::cout << " name = " << name_ << std::endl;
  //std::cout << " objVarEntries = " << objVarEntries.size() << std::endl;
  //std::cout << " prob size = " << prob.size() << std::endl;
  //std::cout << " fixProb size = " << fixProb.size() << std::endl;
  
  if ( !(objVarEntries.size() == prob.size() && objVarEntries.size() == fixProb.size()) ) {
    edm::LogError ("pdfSingleProcessEntryType") << " Number of objVarEntries does not match number of defined probability variables !!";
    assert(0);
  }

  unsigned numVar = objVarEntries.size();
  for ( unsigned iVar = 0; iVar < numVar; ++iVar ) {
    objVarEntryType* objVarEntry = objVarEntries[iVar];
    std::string pdf1dEntryName = std::string(name_).append("_").append(objVarEntry->name_);    
    pdf1dEntryType* pdf1dEntry = new pdf1dEntryType(pdf1dEntryName.data(), objVarEntry, prob[iVar], fixProb[iVar]);
    pdf1dEntries_.push_back(pdf1dEntry);
  }

  std::string normName = std::string(name_).append("_norm");
  std::string normTitle = normName;
  if ( fixNorm ) {
    std::cout << "--> creating RooConstVar with name = " << normName << ":" 
	      << " value = " << norm << "..." << std::endl;
    norm_ = new RooConstVar(normName.data(), normTitle.data(), norm);
  } else {
    std::cout << "--> creating RooRealVar with name = " << normName << ":" 
	      << " start-value = " << norm << "..." << std::endl;
    norm_ = new RooRealVar(normName.data(), normTitle.data(), norm, 0., maxNorm);
  }
  fixNorm_ = fixNorm;

  std::string pdfSingleProcessName = std::string(name_).append("_pdfSingleProcess");
  std::string pdfSingleProcessTitle = pdfSingleProcessName;
  TObjArray pdf_factorCollection;
  for ( std::vector<pdf1dEntryType*>::iterator pdf1dEntry = pdf1dEntries_.begin();
	pdf1dEntry != pdf1dEntries_.end(); ++pdf1dEntry ) {
    pdf_factorCollection.Add((*pdf1dEntry)->pdf1d_);
  }
  std::string pdf_factorName = std::string(name_).append("_args");
  RooArgList pdf_factors(pdf_factorCollection, pdf_factorName.data());
  std::cout << "--> creating RooProdPdf with name = " << pdfSingleProcessName << std::endl;
  pdfSingleProcess_ = new RooProdPdf(pdfSingleProcessName.data(), pdfSingleProcessTitle.data(), pdf_factors);
}

GenMatrixBgEstFit::pdfSingleProcessEntryType::~pdfSingleProcessEntryType()
{
  delete pdfSingleProcess_;

  for ( std::vector<pdf1dEntryType*>::iterator it = pdf1dEntries_.begin();
	it != pdf1dEntries_.end(); ++it ) {
    delete (*it);
  }

  delete norm_;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

GenMatrixBgEstFit::pdfProcessSumEntryType::pdfProcessSumEntryType(const std::string& name, 
								  std::vector<pdfSingleProcessEntryType*>& pdfSingleProcessEntries)
  : name_(name)
{
  //std::cout << "<pdfProcessSumEntryType::pdfProcessSumEntryType>:" << std::endl;
  //std::cout << " name = " << name_ << std::endl;
 
  pdfSingleProcessEntries_ = pdfSingleProcessEntries;

  std::string pdfProcessSumName = std::string(name_).append("_pdfProcessSum");
  std::string pdfProcessSumTitle = pdfProcessSumName;
  TObjArray pdf_summandCollection;
  TObjArray norm_summandCollection;
  for ( std::vector<pdfSingleProcessEntryType*>::iterator pdfSingleProcessEntry = pdfSingleProcessEntries.begin();
	pdfSingleProcessEntry != pdfSingleProcessEntries.end(); ++pdfSingleProcessEntry ) {
    pdf_summandCollection.Add((*pdfSingleProcessEntry)->pdfSingleProcess_);
    norm_summandCollection.Add((*pdfSingleProcessEntry)->norm_);
  }
  std::string pdf_summandName = std::string(name_).append("_args");
  RooArgList pdf_summands(pdf_summandCollection, pdf_summandName.data());
  std::string norm_summandName = std::string(name_).append("_args");
  RooArgList norm_summands(norm_summandCollection, norm_summandName.data());
  std::cout << "--> creating RooAddPdf with name = " << pdfProcessSumName << std::endl;
  pdfProcessSum_ = new RooAddPdf(pdfProcessSumName.data(), pdfProcessSumTitle.data(), pdf_summands, norm_summands); 
}

GenMatrixBgEstFit::pdfProcessSumEntryType::~pdfProcessSumEntryType()
{
  delete pdfProcessSum_;
  for ( std::vector<pdfSingleProcessEntryType*>::iterator it = pdfSingleProcessEntries_.begin();
	it != pdfSingleProcessEntries_.end(); ++it ) {
    delete (*it);
  }
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

GenMatrixBgEstFit::processTreeEntryType::processTreeEntryType(const std::string& processName, 
							      const std::vector<std::string>& fileNames, 
							      const std::string& treeName, const std::string& treeSelection, 
							      std::vector<objVarEntryType*>& objVarEntries,
							      const std::string& branchName_eventWeight)
  : name_(processName)
{
  std::cout << "<processTreeEntryType::processTreeEntryType>:" << std::endl;
  std::cout << " name = " << processName << std::endl;

//--- create chain of all ROOT files 
//    specified in fileNames configuration parameter
  allEventsTree_ = new TChain(treeName.data());
  for ( std::vector<std::string>::const_iterator fileName = fileNames.begin();
	fileName != fileNames.end(); ++fileName ) {
    allEventsTree_->Add(fileName->data());
  }

  std::cout << " entries in Tree = " << allEventsTree_->GetEntries() << std::endl;

//--- copy tree of all events passing selection
  selEventsTree_ = ( treeSelection != "" ) ? allEventsTree_->CopyTree(treeSelection.data()) : allEventsTree_;
  std::cout << " selection = " << treeSelection << std::endl;
  std::cout << " entries passing selection = " << selEventsTree_->GetEntries() << std::endl;

//--- copy tree of all events within phase-space limits
//    considered for background estimation by matrix method
  std::string goodEventsTreeName = std::string(processName).append("_").append("goodEventsTree");
  std::string goodEventsTreeTitle = goodEventsTreeName;
  goodEventsTree_ = new TTree(goodEventsTreeName.data(), goodEventsTreeTitle.data());

  unsigned numVar = objVarEntries.size();
  objVarValues_.resize(numVar);

  for ( unsigned iVar = 0; iVar < numVar; ++iVar ) {
    const std::string& branchName = objVarEntries[iVar]->name_;
    selEventsTree_->SetBranchAddress(branchName.data(), &objVarValues_[iVar]);
    
    std::string leafList = std::string(branchName).append("/F");
    goodEventsTree_->Branch(branchName.data(), &objVarValues_[iVar], leafList.data());
  }

  if ( branchName_eventWeight != "" ) {
    selEventsTree_->SetBranchAddress(branchName_eventWeight.data(), &eventWeight_);
    
    std::string leafList = std::string(branchName_eventWeight).append("/F");
    goodEventsTree_->Branch(branchName_eventWeight.data(), &eventWeight_, leafList.data());
  }
  
  int numEntries = selEventsTree_->GetEntries();
  for ( int iEntry = 0 ; iEntry < numEntries; ++iEntry ) {
    selEventsTree_->GetEvent(iEntry);  
    
    bool isWithinLimits = true;
    for ( unsigned iVar = 0; iVar < numVar; ++iVar ) {
      if ( !(objVarValues_[iVar] >= objVarEntries[iVar]->xMin_ && objVarValues_[iVar] < objVarEntries[iVar]->xMax_) ) 
	isWithinLimits = false;
    }
    
    if ( isWithinLimits ) goodEventsTree_->Fill();
  }

  selEventsTree_->ResetBranchAddresses();

  std::cout << " entries within limits defined by bin-grid = " << goodEventsTree_->GetEntries() << std::endl;

//--- reconfigure tree holding events passing selection and being within phase-space limits 
//    considered for background estimation by matrix method (defined by bin-grid) 
//    for reading
  goodEventsTree_->ResetBranchAddresses();

  for ( unsigned iVar = 0; iVar < numVar; ++iVar ) {
    const std::string& branchName = objVarEntries[iVar]->name_;
    goodEventsTree_->SetBranchAddress(branchName.data(), &objVarValues_[iVar]);
  }

  if ( branchName_eventWeight != "" ) {
    std::cout << "--> reading event weights from branch = " << branchName_eventWeight << std::endl;

    goodEventsTree_->SetBranchAddress(branchName_eventWeight.data(), &eventWeight_);
  }

  std::cout << "done." << std::endl;
}

GenMatrixBgEstFit::processTreeEntryType::~processTreeEntryType()
{
  //std::cout << "<processTreeEntryType::~processTreeEntryType>:" << std::endl;
  //std::cout << " name = " << name_ << std::endl;

  //std::cout << " deleting goodEventsTree..." << std::endl;
  delete goodEventsTree_;
  //std::cout << " deleting selEventsTree..." << std::endl;
  if ( selEventsTree_ != allEventsTree_ ) delete selEventsTree_;
  //std::cout << " deleting allEventsTree..." << std::endl;
  delete allEventsTree_;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

GenMatrixBgEstFit::GenMatrixBgEstFit(const edm::ParameterSet& cfg)
  : cfgError_(0)
{
  //std::cout << "<GenMatrixBgEstFit::GenMatrixBgEstFit>:" << std::endl;

  dataFileNames_ = cfg.getParameter<vstring>("data");

  treeName_ = cfg.getParameter<std::string>("treeName");
  treeSelection_ = ( cfg.exists("treeSelection") ) ? cfg.getParameter<std::string>("treeSelection") : "";

  typedef std::vector<edm::ParameterSet> vParameterSet;
  vParameterSet cfgBranches = cfg.getParameter<edm::ParameterSet>("branches").getParameter<vParameterSet>("config");
  for ( vParameterSet::const_iterator cfgBranch = cfgBranches.begin(); 
	cfgBranch != cfgBranches.end(); ++cfgBranch ) {
    std::string branchName = cfgBranch->getParameter<std::string>("branchName");
    
    branchNames_.push_back(branchName);

    edm::ParameterSet cfgBinning = cfgBranch->getParameter<edm::ParameterSet>("binning");

    typedef std::vector<double> vdouble;
    vdouble binBoundaries = cfgBinning.getParameter<vdouble>("boundaries");
    if ( binBoundaries.size() != 1 ) {
      edm::LogError ("GenMatrixBgEstFit") << " Unsupported number of binBoundaries = " << binBoundaries.size() << " !!";
      cfgError_ = 1;
      continue;
    }

    double xBoundary = binBoundaries[0];
    double xMin = cfgBinning.getParameter<double>("min");
    double xMax = cfgBinning.getParameter<double>("max");
    
    objVarEntryType* objVarEntry = new objVarEntryType(branchName, xBoundary, xMin, xMax);
    objVarEntries_.push_back(objVarEntry);
  }

  branchName_eventWeight_ = "";
  eventWeight_ = 0;
  if ( cfg.exists("branchNameEventWeight") ) {
    branchName_eventWeight_ = cfg.getParameter<std::string>("branchNameEventWeight");

    std::string varEventWeightName = branchName_eventWeight_;
    std::string varEventWeightTitle = varEventWeightName;
    eventWeight_ = new RooRealVar(varEventWeightName.data(), varEventWeightTitle.data(), 1.);
  }

  numVar_ = objVarEntries_.size();
  if ( numVar_ < 2 ) {
    edm::LogError ("GenMatrixBgEstFit") << " Unsupported number of dimensions = " << numVar_ << " !!";
    cfgError_ = 1;
  }

//--- configure RooFit structure
  std::vector<pdfSingleProcessEntryType*> pdfSingleProcessEntries;
  
  edm::ParameterSet cfgProcesses = cfg.getParameter<edm::ParameterSet>("processes");
  vstring processNames = cfgProcesses.getParameterNamesForType<edm::ParameterSet>();
  for ( vstring::const_iterator processName = processNames.begin(); 
	processName != processNames.end(); ++processName ) {
    edm::ParameterSet cfgProcess = cfgProcesses.getParameter<edm::ParameterSet>(*processName);

    vstring fileNames = ( cfgProcess.exists("fileNames") ) ? cfgProcess.getParameter<vstring>("fileNames") : vstring();

//--- estimate probabilities and normalization values for signal/background Monte Carlo processes 
//    from ntuples
    if ( fileNames.size() ) {
      treeData_[*processName] = new processTreeEntryType(*processName, fileNames, treeName_, treeSelection_,
							 objVarEntries_, branchName_eventWeight_);
    }
    
    float norm = 0.;
    if ( cfgProcess.exists("norm") ) {
      cfgProcess.getParameter<double>("norm");
    } else {
      //norm = defaultNorm;
//--- initialize probability values with random number
//    in the range ]sqrt(10),sqrt(10)*1000[      
      norm = TMath::Power(10, 2 + 0.5*getRandom());

      std::cout << (*processName) << ": Norm = " << norm << " assigned randomly." << std::endl;
    }

    bool fixNorm = ( cfgProcess.exists("fixNorm") ) ? cfgProcess.getParameter<bool>("fixNorm") : false;

    //if ( fixNorm && treeData_[*processName] ) {
    if ( treeData_[*processName] ) {
      if ( cfgProcess.exists("norm") ) edm::LogWarning ("GenMatrixBgEstFit") << " Norm constraint twice --> taking value from TTree.";

      norm = compNorm(treeData_[*processName]);

      std::cout << (*processName) << ": Norm = " << norm << " computed from TTree." << std::endl;
    }

    std::vector<float> prob(numVar_);
    std::vector<bool> fixProb(numVar_);
    for ( unsigned iVar = 0; iVar < numVar_; ++iVar ) {

      std::ostringstream cfgProbName;
      cfgProbName << "P" << (iVar + 1);

      prob[iVar] = 0.;
      if ( cfgProcess.exists(cfgProbName.str().data()) ) {
	prob[iVar] = cfgProcess.getParameter<double>(cfgProbName.str().data());
      } else {
	//prob[iVar] = defaultProb;
//--- initialize probability values with random number
//    in the range ]0.2,0.8[
//
//    NOTE: pdf does not depend on observables at all,
//          if default value of 0.5 is used,
//          potentially causing problems with RooFit (?)
//
	prob[iVar] = 0.5 + 0.1*getRandom();

	std::cout << (*processName) << ": " << cfgProbName.str() << " = " << prob[iVar] << " assigned randomly." << std::endl;
      }

      std::string cfgFixProbName = std::string("fix").append(cfgProbName.str());
      fixProb[iVar] = ( cfgProcess.exists(cfgFixProbName.data()) ) ? 
	cfgProcess.getParameter<bool>(cfgFixProbName.data()) : false;

      //if ( fixProb[iVar] && treeData_[*processName] ) {
      if ( treeData_[*processName] ) {
	if ( cfgProcess.exists(cfgProbName.str().data()) ) 
	  edm::LogWarning ("GenMatrixBgEstFit") << " Probability " << cfgProbName.str() << " constraint twice"
						<< " --> taking value from TTree.";

	prob[iVar] = compProb(treeData_[*processName], iVar);

	std::cout << (*processName) << ": " << cfgProbName.str() << " = " << prob[iVar] << " computed from TTree," 
		  << " branch name = " << objVarEntries_[iVar]->name_ << "." << std::endl;
      }
    }

    pdfSingleProcessEntryType* pdfSingleProcessEntry = new pdfSingleProcessEntryType(*processName, objVarEntries_, 
										     prob, fixProb, norm, fixNorm);

    edm::ParameterSet cfgDrawOptions_process = cfgProcess.getParameter<edm::ParameterSet>("drawOptions");
    pdfSingleProcessEntry->lineColor_ = cfgDrawOptions_process.getParameter<int>("lineColor");
    pdfSingleProcessEntry->lineStyle_ = cfgDrawOptions_process.getParameter<int>("lineStyle");
    pdfSingleProcessEntry->lineWidth_ = cfgDrawOptions_process.getParameter<int>("lineWidth");

    pdfSingleProcessEntries.push_back(pdfSingleProcessEntry);
  }

  dataSet_ = 0;

  std::string pdfProcessSumName = "pdfProcessSum";
  pdfProcessSumEntry_ = new pdfProcessSumEntryType(pdfProcessSumName, pdfSingleProcessEntries);

  edm::ParameterSet cfgBinGrid = cfg.getParameter<edm::ParameterSet>("branches");
  binGrid_ = new BinGrid(cfgBinGrid); 

  edm::ParameterSet cfgControlPlots = cfg.getParameter<edm::ParameterSet>("output").getParameter<edm::ParameterSet>("controlPlots");
  controlPlotsFileName_ = cfgControlPlots.getParameter<std::string>("fileName");

  edm::ParameterSet cfgScaleFactors = cfg.getParameter<edm::ParameterSet>("output").getParameter<edm::ParameterSet>("scaleFactors");
  scaleFactorFileName_ = cfgScaleFactors.getParameter<std::string>("fileName");
  scaleFactorSignalRegion_ = cfgScaleFactors.getParameter<vdouble>("signalRegion");
}

GenMatrixBgEstFit::~GenMatrixBgEstFit()
{
  std::cout << "<GenMatrixBgEstFit::~GenMatrixBgEstFit>:" << std::endl;

//--- delete RooFit objects
  //std::cout << " deleting dataSet..." << std::endl;
  delete dataSet_;

  //std::cout << " deleting pdfProcessSumEntry..." << std::endl;
  delete pdfProcessSumEntry_;

  //std::cout << " deleting objVarEntries..." << std::endl;
  for ( std::vector<objVarEntryType*>::iterator it = objVarEntries_.begin();
	it != objVarEntries_.end(); ++it ) {
    delete (*it);
  }

//--- close files containing ntuples 
//    for signal/background Monte Carlo processes and (pseudo)data
  for ( std::map<std::string, processTreeEntryType*>::iterator it = treeData_.begin();
	it != treeData_.end(); ++it ) {
    delete it->second;
  }

  delete binGrid_;
}

void GenMatrixBgEstFit::endJob()
{
  //std::cout << "<GenMatrixBgEstFit::endJob>:" << std::endl; 

//--- check that configuration parameters contain no errors
  if ( cfgError_ ) {
    edm::LogError ("GenMatrixBgEstFit::endJob") << " Error in Configuration ParameterSet --> skipping !!";
    return;
  }

//--- extract tree entries that are within the [xMin, xMax[ limits
//    specified for each branch
  treeData_["data"] = new processTreeEntryType("data", dataFileNames_, treeName_, treeSelection_,
					       objVarEntries_, branchName_eventWeight_);

//--- create RooFit dataset;
//    fit pdf
  std::string dataSetName = "dataSet";
  std::string dataSetTitle = dataSetName;
  TObjArray dataSet_variableCollection;
  for ( std::vector<objVarEntryType*>::const_iterator objVarEntry = objVarEntries_.begin();
	objVarEntry != objVarEntries_.end(); ++objVarEntry ) {
    dataSet_variableCollection.Add((*objVarEntry)->x_);
  }
  if ( eventWeight_ ) dataSet_variableCollection.Add(eventWeight_);
  std::string dataSet_variableName = "dataSet_args";
  RooArgSet dataSet_variables(dataSet_variableCollection, dataSet_variableName.data());
  dataSet_ = new RooDataSet(dataSetName.data(), dataSetTitle.data(), treeData_["data"]->goodEventsTree_, dataSet_variables);
  //dataSet_ = new RooDataSet(dataSetName.data(), dataSetTitle.data(), treeData_["data"]->allEventsTree_, dataSet_variables);
  if ( eventWeight_ ) dataSet_->setWeightVar(*eventWeight_);

  RooAbsPdf* model = pdfProcessSumEntry_->pdfProcessSum_;

  std::cout << ">>> RootFit model used for generalized Matrix method Fit <<<" << std::endl;
  model->printCompactTree();

  std::cout << ">>> RootFit Variables <<<" << std::endl;
  model->getVariables()->Print("v");

  std::cout << ">>> RootFit Parameters <<<" << std::endl;
  model->getParameters(dataSet_)->Print("v");

  std::cout << ">>> RootFit Observables <<<" << std::endl;
  model->getObservables(dataSet_)->Print("v");

  model->fitTo(*dataSet_, RooFit::Extended());
  //model->fitTo(*dataSet_);

//--- print-out fit results
  print(std::cout);

//--- produce plots of different signal and background processes
//    using scale factors determined by fit
//    compared to distribution of (pseudo)data
  if ( controlPlotsFileName_ != "" ) {
    for ( unsigned iVar = 0; iVar < numVar_; ++iVar ) {
      const objVarEntryType* objVarEntry = objVarEntries_[iVar];

      std::ostringstream variableName;
      variableName << "variable" << " " << (iVar + 1);
      std::ostringstream variableTitle;
      variableTitle << "var_{" << (iVar + 1) << "}";

      if ( controlPlotsFileName_.find(".") == std::string::npos ||
	   controlPlotsFileName_.find(".") == (controlPlotsFileName_.length() - 1) ) {
	edm::LogError ("GenMatrixBgEstFit::endJob") << " Invalid format for fileName = " << controlPlotsFileName_
						    << " --> skipping !!";
	return;
      } 

      size_t posSeparator = controlPlotsFileName_.find_last_of(".");
      std::ostringstream outputFileName;
      outputFileName << std::string(controlPlotsFileName_, 0, posSeparator)
		     << "_var" << (iVar + 1)
		     << std::string(controlPlotsFileName_, posSeparator);
      std::cout << " outputFileName = " << outputFileName.str() << std::endl;

      makeControlPlot(objVarEntry->x_, variableName.str(), variableTitle.str(), outputFileName.str());
    }
  }

//--- produce python configuration file 
//    containing scale factors/normalization constants determined by fit
  if ( scaleFactorFileName_ != "" ) makeScaleFactors();
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

void GenMatrixBgEstFit::print(std::ostream& stream)
{
  stream << "<GenMatrixBgEstFit::print>:" << std::endl;

  unsigned numVar = objVarEntries_.size();

  stream << "Fit Parameter:" << std::endl;
  std::vector<pdfSingleProcessEntryType*>& pdfSingleProcessEntries = pdfProcessSumEntry_->pdfSingleProcessEntries_;
  for ( std::vector<pdfSingleProcessEntryType*>::const_iterator pdfSingleProcessEntry = pdfSingleProcessEntries.begin(); 
	pdfSingleProcessEntry != pdfSingleProcessEntries.end(); ++pdfSingleProcessEntry ) {
    stream << " Process = " << (*pdfSingleProcessEntry)->name_ << std::endl;
    std::string fixNorm_string = ( (*pdfSingleProcessEntry)->fixNorm_ ) ? "(fixed)" : "";
    stream << "  norm = " << (*pdfSingleProcessEntry)->norm_->getVal() << " " << fixNorm_string << std::endl;
    std::vector<pdf1dEntryType*>& pdf1dEntries = (*pdfSingleProcessEntry)->pdf1dEntries_;
    for ( unsigned iVar = 0; iVar < numVar; ++iVar ) {
      std::ostringstream probName;
      probName << "P" << (iVar + 1);
      std::string fixProb_string = ( pdf1dEntries[iVar]->fixProb_ ) ? "(fixed)" : "";
      const Float_t* binBoundaries = pdf1dEntries[iVar]->pdf1d_->getBinBoundaries();
      double prob1, prob2;
      getProbability(pdf1dEntries[iVar]->prob1_->getVal(), pdf1dEntries[iVar]->prob2_->getVal(), binBoundaries, prob1, prob2);
      //stream << "  " << probName.str() << " = " << pdf1dEntries[iVar]->prob1_->getVal() << " " << fixProb_string << std::endl;
      stream << "  " << probName.str() << " = " << prob1 << " " << fixProb_string << std::endl;
     }
  }

  stream << "Contributions to Regions:" << std::endl;
  unsigned numRegions = TMath::Nint(TMath::Power(2, numVar));
  for ( unsigned iRegion = 0; iRegion < numRegions; ++iRegion ) {
    std::ostringstream regionDescription1;
    std::ostringstream regionDescription2;
    regionDescription2 << "Norm * ";
    for ( unsigned iVar = 0; iVar < numVar; ++iVar ) {
      unsigned bitValue = TMath::Nint(TMath::Power(2, iVar));

      std::string operator_string;
      std::ostringstream probDescription;
      if ( iRegion & bitValue ) {
	operator_string = "<";
	probDescription << "(1 - P" << (iVar + 1) << ")";
      } else {
	operator_string = ">=";
	probDescription << "P" << (iVar + 1);
      }

      regionDescription1 << objVarEntries_[iVar]->name_ << " " << operator_string << " " << objVarEntries_[iVar]->xBoundary_;
      regionDescription2 << probDescription.str();
      if ( iVar < (numVar - 1) ) {
	regionDescription1 << " && "; 
	regionDescription2 << " * "; 
      }
    }

    stream << " " << regionDescription1.str() << std::endl;
    stream << " " << regionDescription2.str() << std::endl;

    std::vector<pdfSingleProcessEntryType*>& pdfSingleProcessEntries = pdfProcessSumEntry_->pdfSingleProcessEntries_;
    for ( std::vector<pdfSingleProcessEntryType*>::const_iterator pdfSingleProcessEntry = pdfSingleProcessEntries.begin(); 
	  pdfSingleProcessEntry != pdfSingleProcessEntries.end(); ++pdfSingleProcessEntry ) {
      float processContribution = (*pdfSingleProcessEntry)->norm_->getVal();
      std::vector<pdf1dEntryType*>& pdf1dEntries = (*pdfSingleProcessEntry)->pdf1dEntries_;
      for ( unsigned iVar = 0; iVar < numVar; ++iVar ) {
	unsigned bitValue = TMath::Nint(TMath::Power(2, iVar));

	const Float_t* binBoundaries = pdf1dEntries[iVar]->pdf1d_->getBinBoundaries();
	double prob1, prob2;
	getProbability(pdf1dEntries[iVar]->prob1_->getVal(), pdf1dEntries[iVar]->prob2_->getVal(), binBoundaries, prob1, prob2);

	if ( iRegion & bitValue ) 
	  //processContribution *= (*pdfSingleProcessEntry)->pdf1dEntries_[iVar]->prob1_->getVal();
	  processContribution *= prob1;
	else 
	  //processContribution *= (1 - (*pdfSingleProcessEntry)->pdf1dEntries_[iVar]->prob1_->getVal());
	  processContribution *= prob2;
      }

      stream << "  " << (*pdfSingleProcessEntry)->name_ << " = " 
	     << std::setprecision(3) << std::fixed << processContribution << std::endl;
    }
  }
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

void GenMatrixBgEstFit::makeControlPlot(const RooRealVar* x, 
					const std::string& variableName, const std::string& variableTitle,
					const std::string& outputFileName)
{
//--- stop ROOT from opening X-window for canvas output
//    (in order to be able to run in batch mode) 
  gROOT->SetBatch(true);

  TCanvas canvas("GenMatrixBgEstFit", "GenMatrixBgEstFit", defaultCanvasSizeX, defaultCanvasSizeY);
  canvas.SetFillColor(10);

  RooPlot* plotFrame = x->frame();
  std::string plotTitle = std::string("GenMatrixBgEstFit - Control Distribution of ").append(variableName);
  plotFrame->SetTitle(plotTitle.data());
  dataSet_->plotOn(plotFrame, RooFit::MarkerColor(kBlack), RooFit::MarkerStyle(2));
  const std::vector<pdfSingleProcessEntryType*>& pdfSingleProcessEntries = pdfProcessSumEntry_->pdfSingleProcessEntries_;
  for ( std::vector<pdfSingleProcessEntryType*>::const_iterator pdfSingleProcessEntry = pdfSingleProcessEntries.begin();
	pdfSingleProcessEntry != pdfSingleProcessEntries.end(); ++pdfSingleProcessEntry ) {
    std::string componentName = (*pdfSingleProcessEntry)->pdfSingleProcess_->GetName();
    pdfProcessSumEntry_->pdfProcessSum_->plotOn(plotFrame, RooFit::Components(componentName.data()), 
						RooFit::LineColor((*pdfSingleProcessEntry)->lineColor_), 
						RooFit::LineStyle((*pdfSingleProcessEntry)->lineStyle_), 
						RooFit::LineWidth((*pdfSingleProcessEntry)->lineWidth_));
  }
  pdfProcessSumEntry_->pdfProcessSum_->plotOn(plotFrame, RooFit::LineColor(kBlack), RooFit::LineStyle(kSolid), RooFit::LineWidth(2));
  plotFrame->SetXTitle(variableTitle.data());
  //plotFrame->SetTitleOffset(1.2, "X");
  plotFrame->SetYTitle("");
  //plotFrame->SetTitleOffset(1.2, "Y");
  plotFrame->Draw();
  
  canvas.Update();
  
  canvas.Print(outputFileName.data());
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

void GenMatrixBgEstFit::makeScaleFactors()
{
//--- check that dimensionality of point specifying signal region
//    matches number of observables
  if ( scaleFactorSignalRegion_.size() != numVar_ ) {
    edm::LogError ("makeScaleFactors") << " Dimensionality of point specifying signal region = " << scaleFactorSignalRegion_.size()
				       << " does not match Number of observables = " << numVar_ << " --> skipping !!";
    return;
  }

//--- check that point specifying signal point
//    is within defined bin-grid
  unsigned numRegions = TMath::Nint(TMath::Power(2, numVar_));
  unsigned iRegion = binGrid_->binNumber(scaleFactorSignalRegion_);
  if ( !(iRegion >= 1 && iRegion <= numRegions) ) {
    edm::LogError ("makeScaleFactors") << " Point specifying signal region not within defined bin-grid --> skipping !!";
    return;
  }

//--- open output file
  std::ostream* outputFile = new std::ofstream(scaleFactorFileName_.data(), std::ios::out);

  std::vector<pdfSingleProcessEntryType*>& pdfSingleProcessEntries = pdfProcessSumEntry_->pdfSingleProcessEntries_;
  for ( std::vector<pdfSingleProcessEntryType*>::const_iterator pdfSingleProcessEntry = pdfSingleProcessEntries.begin(); 
	pdfSingleProcessEntry != pdfSingleProcessEntries.end(); ++pdfSingleProcessEntry ) {

//--- compute normalization constant
    float processContribution = (*pdfSingleProcessEntry)->norm_->getVal();
    std::vector<pdf1dEntryType*>& pdf1dEntries = (*pdfSingleProcessEntry)->pdf1dEntries_;
    for ( unsigned iVar = 0; iVar < numVar_; ++iVar ) {
      unsigned bitValue = TMath::Nint(TMath::Power(2, iVar));

      const Float_t* binBoundaries = pdf1dEntries[iVar]->pdf1d_->getBinBoundaries();
      double prob1, prob2;
      getProbability(pdf1dEntries[iVar]->prob1_->getVal(), pdf1dEntries[iVar]->prob2_->getVal(), binBoundaries, prob1, prob2);

      if ( iRegion & bitValue ) 
	processContribution *= prob1;
      else 
	processContribution *= prob2;
    }

//--- write normalization constant for process into output file
    (*outputFile) << (*pdfSingleProcessEntry)->name_ << ".normalization = cms.double(" 
		  << std::setprecision(3) << std::fixed << processContribution << ")" << std::endl;
  }

//--- close output file
  delete outputFile;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

double GenMatrixBgEstFit::compProb(processTreeEntryType* processTreeEntry, unsigned iVar)
{
  //std::cout << "<GenMatrixBgEstFit::compProb>:" << std::endl;
  //std::cout << " goodEventsTreeName = " << processTreeEntry->goodEventsTree_->GetName() << std::endl;
  //std::cout << " xBoundary = " << objVarEntries_[iVar]->xBoundary_ << std::endl;
  //std::cout << " xMin = " << objVarEntries_[iVar]->xMin_ << std::endl;
  //std::cout << " xMax = " << objVarEntries_[iVar]->xMax_ << std::endl;

  double numValues_belowXboundary = 0.;
  double numValues_aboveXboundary = 0.;

  int numEntries = processTreeEntry->goodEventsTree_->GetEntries();
  for ( int iEntry = 0; iEntry < numEntries; ++iEntry ) {
    processTreeEntry->goodEventsTree_->GetEvent(iEntry);  

    if ( processTreeEntry->objVarValues_[iVar] >= objVarEntries_[iVar]->xMin_ && 
	 processTreeEntry->objVarValues_[iVar] <  objVarEntries_[iVar]->xBoundary_ ) 
      numValues_belowXboundary += processTreeEntry->eventWeight_;
    else if ( processTreeEntry->objVarValues_[iVar] >= objVarEntries_[iVar]->xBoundary_ &&
	      processTreeEntry->objVarValues_[iVar] <  objVarEntries_[iVar]->xMax_ ) 
      numValues_aboveXboundary += processTreeEntry->eventWeight_;
  }

  return numValues_belowXboundary/(numValues_belowXboundary + numValues_aboveXboundary);
}

double GenMatrixBgEstFit::compNorm(processTreeEntryType* processTreeEntry)
{
  //std::cout << "<GenMatrixBgEstFit::compNorm>:" << std::endl;
  //std::cout << " goodEventsTreeName = " << processTreeEntry->goodEventsTree_->GetName() << std::endl;
  //std::cout << " numEntries = " << processTreeEntry->goodEventsTree_->GetEntries() << std::endl;

  if ( branchName_eventWeight_ == "" ) return processTreeEntry->goodEventsTree_->GetEntries();

  double norm = 0.;

  int numEntries = processTreeEntry->goodEventsTree_->GetEntries();
  for ( int iEntry = 0; iEntry < numEntries; ++iEntry ) {
    processTreeEntry->goodEventsTree_->GetEvent(iEntry);

    norm += processTreeEntry->eventWeight_;
  }

  return norm;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(GenMatrixBgEstFit);
