#include "TauAnalysis/BgEstimationTools/plugins/TF1Wrapper.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <string>

TF1Wrapper::TF1Wrapper(const edm::ParameterSet& cfg)
  : TF1WrapperBase(cfg) 
{
  std::string pluginName = ( cfg.exists("pluginName") ) ? cfg.getParameter<std::string>("pluginName") : "TF1Wrapper::tf1_";
  
  std::string formula = cfg.getParameter<std::string>("formula");
  
  double xMin = cfg.getParameter<double>("xMin");
  double xMax = cfg.getParameter<double>("xMax");
  
  tf1_ = new TF1(pluginName.data(), formula.data(), xMin, xMax);

  setTF1Parameter(cfg.getParameter<edm::ParameterSet>("parameter"));
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(TF1WrapperPluginFactory, TF1Wrapper, "TF1Wrapper");
