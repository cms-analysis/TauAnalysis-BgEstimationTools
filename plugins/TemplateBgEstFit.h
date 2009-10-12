#ifndef TauAnalysis_BgEstimationTools_TemplateBgEstFit_h  
#define TauAnalysis_BgEstimationTools_TemplateBgEstFit_h

/** \class TemplateBgEstFit
 *
 * Estimate contribution of signal and background processes
 * to final event sample by fitting shape "templates" for different processes
 * to distribution observed in data
 * (class implements "template" method for data-driven background estimation)
 *
 * NOTE: The TemplateBgEstFit class is capable of fitting distributions of up to three different observables simultaneously,
 *       determining the contribution of signal and background processes as the normalization parameters of the templates
 *       that best fit the combination of all observed distributions.
 *       
 *       In case of more than one distribution fitted, the fit is implemented by constructing an auxiliary product PDF,
 *       obtained by multiplying the different one-dimensional distributions to be fitted.
 *
 *       This multiplication is done for purely technical reasons only 
 *       (in order to implement the solution to the problem 
 *        "fit one, two or three one-dimensional distributions by a sum of template PDFs
 *         the normalization parameters of which are to be determined by the fit in such a way
 *         that for each signal/background process the same normalization parameter fits all distributions"
 *        in a way that is supported by the RooFit package used to implement the technical aspects of the fit)
 *       and does in particular NOT mean that any assumption is made 
 *       that the observables used in the fit are uncorrelated/independent of each other !!
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.2.2.2 $
 *
 * $Id: TemplateBgEstFit.h,v 1.2.2.2 2009/09/16 11:47:00 veelken Exp $
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

#include <string>

class TemplateBgEstFit : public edm::EDAnalyzer
{
  struct dataDistr1dType
  {
    dataDistr1dType(const std::string&, const std::string&, const std::string&, RooRealVar*, bool);
    virtual ~dataDistr1dType();
    virtual void initialize();
    void buildFitData();
    virtual void fluctuate(bool, bool);
    std::string processName_;
    std::string varName_;
    std::string meName_;
    MonitorElement* me_;
    RooRealVar* xRef_;
    bool cutUnfittedRegion_;
    TH1* histogram_;    
    std::string dataHistName_;
    RooDataHist* dataHist_;
    TH1* fluctHistogram_;
    int error_;
  };

  struct dataDistrNdType 
  {
    dataDistrNdType(int, bool);
    ~dataDistrNdType();
    void addElement(const std::string&, RooRealVar*, const std::string&);
    void initialize();
    void buildFitData();
    void fluctuate(bool, bool);
    unsigned numDimensions_;
    std::vector<std::string> varNames_;
    std::map<std::string, dataDistr1dType*> dataEntries1d_;
    int fitMode_;
    bool cutUnfittedRegion_;
    RooDataHist* fitData_;
    double normCorrFactor_;
    TH1* auxHistogram_;
    int error_;
  };

  struct sysFluctDefType
  {
    std::string fluctName_;
    std::string meName_;
    MonitorElement* me_;
    double pullRMS_;
    double pullMin_;
    double pullMax_;
    int fluctMode_;
  };

  struct modelTemplate1dType : public dataDistr1dType
  {
    modelTemplate1dType(const std::string&, const std::string&, const std::string&, RooRealVar*, bool, bool, const edm::ParameterSet&);
    virtual ~modelTemplate1dType();
    void initialize();
    void buildPdf();
    void fluctuate(bool, bool);
    std::string pdf1dName_;
    RooAbsPdf* pdf1d_;
    bool applySmoothing_;
    edm::ParameterSet cfgSmoothing_;
    TF1WrapperBase* auxTF1Wrapper_;
    std::vector<sysFluctDefType> sysErrFluctuations_;
  };

  struct modelTemplateNdType
  {
    modelTemplateNdType(const std::string&, const edm::ParameterSet&, int, bool, bool, double, double);
    ~modelTemplateNdType();
    void addElement(const std::string&, RooRealVar*, const std::string&, bool, const edm::ParameterSet&);
    void initialize();
    void buildPdf();
    void fluctuate(bool, bool);
    std::string processName_;
    bool applyNormConstraint_;
    RooAbsPdf* pdfNormConstraint_;
    RooConstVar* meanNormConstraint_;
    RooConstVar* sigmaNormConstraint_;
    std::string pdfName_;
    RooAbsPdf* pdf_;
    bool pdfIsOwned_;
    RooRealVar* norm_;
    unsigned numDimensions_;
    std::vector<std::string> varNames_;
    std::map<std::string, modelTemplate1dType*> processEntries1d_;
    int fitMode_;
    bool cutUnfittedRegion_;
    double normCorrFactor_;
    int lineColor_;
    int lineStyle_;
    int lineWidth_;
    TH1* auxHistogram_;
    RooDataHist* auxDataHist_;
    RooAbsPdf* auxPdf_;
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
  std::vector<std::string> processNames_;
  std::vector<std::string> varNames_;

  typedef std::map<std::string, modelTemplateNdType*> processEntryMap;
  processEntryMap processEntries_;

  dataDistrNdType* dataEntry_;

  typedef std::map<std::string, RooRealVar*> realVarMap;
  realVarMap x_;

  RooAbsPdf* fitModel_;
  RooFitResult* fitResult_;
  TVectorD fitResultMean_;
  TMatrixD fitResultCov_;

  int error_;
};

#endif  


