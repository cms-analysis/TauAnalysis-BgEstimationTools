#ifndef TauAnalysis_BgEstimationTools_GenMatrix1dPdf_h
#define TauAnalysis_BgEstimationTools_GenMatrix1dPdf_h

/** \class GenMatrix1dPdf
 *
 * Auxiliary class for RooFit based background estimation,
 * representing probability for observable X_i to be in "signal-like" region;
 * derived from RooAbsPdf
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: GenMatrix1dPdf.h,v 1.1 2009/06/11 07:23:28 veelken Exp $
 *
 */

#include <TObjArray.h>

#include <RooAbsPdf.h>
#include <RooRealVar.h>
#include <RooRealProxy.h>

class GenMatrix1dPdf : public RooAbsPdf 
{
 public:
  GenMatrix1dPdf(const char*, const char*, RooAbsReal&, Int_t numBins, Double_t* binBoundaries, TCollection&);
  GenMatrix1dPdf(const GenMatrix1dPdf&, const char* = 0);
  virtual TObject* clone(const char* newName) const { return new GenMatrix1dPdf(*this, newName); }
  virtual ~GenMatrix1dPdf();

  int getNumBins() const { return binProbabilities_.GetEntries(); }
  const Double_t* getBinBoundaries() const { return binBoundaries_; }

 protected:
  RooRealProxy x_;

  /// "left" bin-edges (following the convention of ROOT's TH1 class); 
  /// size of array equals number of bins + 1
  Double_t* binBoundaries_; 

  /// (average) probability value for each bin,
  /// stored as type RooAbsArg*
  TObjArray binProbabilities_;      

  Double_t evaluate() const;
};

#endif
