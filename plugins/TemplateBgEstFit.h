#ifndef TauAnalysis_BgEstimationTools_TemplateBgEstFit_h  
#define TauAnalysis_BgEstimationTools_TemplateBgEstFit_h

/** \class TemplateBgEstFit
 *
 * Estimate contribution of signal and background processes to final event sample 
 * by fitting shape "templates" for different processes to distribution observed in data
 * (class implements "template" method for data-driven background estimation)
 *
 * NOTE: The TemplateBgEstFit class is capable of fitting distributions of multiple different observables simultaneously,
 *       determining the contribution of signal and background processes as the normalization parameters of the templates
 *       that best fit the combination of all observed distributions.
 *
 *       A "workaround" is used to handle statistical uncertainties on the templates,
 *       which are not directly supported by the likelihood function build by the RooFit package
 *       that is used by the TemplateBgEstFit class internally to perform fits.
 *
 *       In order to avoid systematically underestimating (overestimating) template shapes with large (small)
 *       statistical uncertainties (cf. http://root.cern.ch/phpBB2/viewtopic.php?p=38815#38815 )
 *       the template shapes determined from background enriched (control) data samples
 *       can by either parametrized by a "smooth" function prior to the fit,
 *       or by fitting the distribution observed in the final event sample simultaneously with the templates
 *       (cf. https://hypernews.cern.ch/HyperNews/CMS/get/ewk/278.html for more details about and discussion of this technique,
 *        originating from BaBar physics analyses)
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.11 $
 *
 * $Id: TemplateBgEstFit.h,v 1.11 2009/10/25 15:47:13 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQMServices/Core/interface/MonitorElement.h"

#include "TauAnalysis/BgEstimationTools/interface/TF1WrapperBase.h"

#include <TVectorD.h>
#include <TMatrixD.h>

#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooRealVar.h>
#include <RooConstVar.h>
#include <RooAbsPdf.h>
#include <RooFitResult.h>
#include <RooCategory.h>
#include <RooSimultaneous.h>

#include <vector>
#include <string>
#include <map>

class TemplateBgEstFit : public edm::EDAnalyzer
{
  struct data1dType
  {
    data1dType(const std::string&, const std::string&, const std::string&, RooRealVar*, bool);
    virtual ~data1dType();
    virtual void initialize();
    virtual void fluctuate(bool, bool);
    std::string processName_;
    std::string varName_;
    std::string meName_;
    MonitorElement* me_;
    RooRealVar* xRef_;
    bool cutUnfittedRegion_;
    TH1* histogram_;    
    TH1* fluctHistogram_;
    int error_;
  };

  struct dataNdType 
  {
    dataNdType(bool);
    ~dataNdType();
    void addVar(const std::string&, RooRealVar*, const std::string&);
    void initialize();
    void fluctuate(bool, bool);
    unsigned numDimensions_;
    std::vector<std::string> varNames_;
    std::map<std::string, data1dType*> data1dEntries_;
    bool cutUnfittedRegion_;
    double normCorrFactor_;
    int error_;
  };

  struct model1dType : public data1dType
  {
    struct sysErrFluctType
    {
      std::string fluctName_;
      std::string meName_;
      MonitorElement* me_;
      double pullRMS_;
      double pullMin_;
      double pullMax_;
      int fluctMode_;
    };

    model1dType(const std::string&, const std::string&, const std::string&, RooRealVar*, bool, bool, const edm::ParameterSet&);
    virtual ~model1dType();
    void initialize();
    void fluctuate(bool, bool);
    void buildPdf();
    std::string pdfName_;
    RooAbsPdf* pdf_;
    bool applySmoothing_;
    edm::ParameterSet cfgSmoothing_;
    TF1WrapperBase* auxTF1Wrapper_;
    TArrayD* pdfBinning_;
    TObjArray* pdfCoeffCollection_;
    RooArgList* pdfCoeffArgs_;
    std::vector<sysErrFluctType> sysErrFluctuations_;
  };

  struct modelNdType
  {
    modelNdType(const std::string&, const edm::ParameterSet&, bool, bool, double, double);
    ~modelNdType();
    void addVar(const std::string&, RooRealVar*, const std::string&, bool, const edm::ParameterSet&);
    void initialize();
    void fluctuate(bool, bool);   
    std::string processName_;
    bool applyNormConstraint_;
    RooAbsPdf* pdfNormConstraint_;
    RooConstVar* meanNormConstraint_;
    RooConstVar* sigmaNormConstraint_;
    RooRealVar* norm_;
    bool cutUnfittedRegion_;
    double normCorrFactor_;
    unsigned numDimensions_;
    std::vector<std::string> varNames_;
    std::map<std::string, model1dType*> model1dEntries_;
    int lineColor_;
    int lineStyle_;
    int lineWidth_;
    int error_;
  };

 public:
  
  explicit TemplateBgEstFit(const edm::ParameterSet&);
  ~TemplateBgEstFit();
  
 private:

  void beginJob(const edm::EventSetup&) {}
  void analyze(const edm::Event&, const edm::EventSetup&) {}
  void endJob();

//--- private auxiliary functions
  void buildFitData();
  void buildFitModel();

  void fit(bool, int, int);

  void print(std::ostream& stream);

  void makeControlPlots();
  void makeControlPlotsSmoothing();
  typedef std::vector<std::string> vstring;
  void makeControlPlotsCovariance(TVectorD, TVectorD, TMatrixD, const vstring&, const std::string&, const char*);
  void makeControlPlotsObsDistribution();

  double compChi2red();

  void estimateUncertainties(bool, bool, int, double, const char*, int, bool);

//--- configuration parameters
  int fitMode_;
  
  bool cutUnfittedRegion_;

  int printLevel_;
  bool printWarnings_;

  std::string controlPlotsFileName_;

  int statErrNumSamplings_;
  double statErrChi2redMax_;
  int statErrPrintLevel_;
  bool statErrPrintWarnings_;

  int sysErrNumSamplings_;
  double sysErrChi2redMax_;
  int sysErrPrintLevel_;
  bool sysErrPrintWarnings_;
  
//--- internal data-members for handling histograms
//    and performing fit
  typedef std::vector<std::string> vstring;
  vstring processNames_;
  vstring varNames_;
  
  dataNdType* dataEntry_;
  RooAbsData* fitData_;

  typedef std::map<std::string, modelNdType*> modelEntryMap;
  modelEntryMap modelEntries_;
  typedef std::map<std::string, RooAbsPdf*> pdfMap;
  pdfMap pdfModelSums_;
  typedef std::map<std::string, RooAbsReal*> normMap;
  normMap normTemplateShapes_;
  pdfMap pdfTemplateShapeSums_;

  //RooSimultaneous* fitModel_;
  RooAbsPdf* fitModel_;

  RooCategory* fitCategories_;

  typedef std::map<std::string, RooRealVar*> realVarMap;
  realVarMap x_;

  RooFitResult* fitResult_;
  TVectorD fitResultMean_;
  TMatrixD fitResultCov_;

  int error_;
};

#endif  


