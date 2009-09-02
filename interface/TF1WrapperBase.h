#ifndef TauAnalysis_BgEstimationTools_TF1WrapperBase_h
#define TauAnalysis_BgEstimationTools_TF1WrapperBase_h

/** \class TF1WrapperBase
 *
 * Base-class for creating TF1 objects via plugin factory
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: TF1WrapperBase.h,v 1.1 2009/09/01 14:17:42 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <TF1.h>

#include <map>

class TF1WrapperBase
{
 public:
  // constructor 
  explicit TF1WrapperBase(const edm::ParameterSet&);
  
  // destructor
  virtual ~TF1WrapperBase();

  // set parameter values of TF1 object to initial values
  void reinitializeTF1Parameter();

  // method for accessing TF1 object
  const TF1* getTF1() const { return ( !cfgError_ ) ? tf1_ : 0; }
  TF1* getTF1() { return ( !cfgError_ ) ? tf1_ : 0; }
  
 protected:
  // auxiliary function to create TF1 object
  // (enforce function to be implemented in derrived classes,
  //  cannot make it a pure virtual function, because it is called from the TF1WrapperBase constructor,
  //  so use an assert instead)
  virtual void makeTF1(const edm::ParameterSet&) { assert(0); } 

  // auxiliary function to initialize parameter of TF1 object
  virtual void setTF1Parameter(const edm::ParameterSet&);

  // set of initial parameter values of TF1 object
  // (needed by reinitializeTF1Parameter function)
  std::map<int, double> tf1InitialValues_;

  // pointer to TF1 object
  TF1* tf1_;

  int cfgError_;
};

#include "FWCore/PluginManager/interface/PluginFactory.h"

typedef edmplugin::PluginFactory<TF1WrapperBase* (const edm::ParameterSet&)> TF1WrapperPluginFactory;

#endif  

