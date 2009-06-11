#ifndef TauAnalysis_BgEstimationTools_GenMatrixBgEst1dPdf_h
#define TauAnalysis_BgEstimationTools_GenMatrixBgEst1dPdf_h

/** \class GenMatrixBgEst1dPdf
 *
 * Auxiliary class for RooFit based background estimation,
 * representing probability for observable X_i to be in "signal-like" region;
 * derived from RooAbsPdf
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: GenMatrixBgEst1dPdf.h,v 1.1 2009/02/04 15:53:56 veelken Exp $
 *
 */

#include <TObjArray.h>

#include <RooAbsPdf.h>
#include <RooRealVar.h>
#include <RooRealProxy.h>

class GenMatrixBgEst1dPdf : public RooAbsPdf 
{
 public:
  GenMatrixBgEst1dPdf(const char*, const char*, RooAbsReal&, Int_t numBins, Float_t* binBoundaries, TCollection&);
  GenMatrixBgEst1dPdf(const GenMatrixBgEst1dPdf&, const char* = 0);
  virtual TObject* clone(const char* newName) const { return new GenMatrixBgEst1dPdf(*this, newName); }
  virtual ~GenMatrixBgEst1dPdf();

  int getNumBins() const { return binProbabilities_.GetEntries(); }
  const Float_t* getBinBoundaries() const { return binBoundaries_; }

 protected:
  RooRealProxy x_;

  /// "left" bin-edges (following the convention of ROOT's TH1 class); 
  /// size of array equals number of bins + 1
  Float_t* binBoundaries_; 

  /// (average) probability value for each bin,
  /// stored as type RooAbsArg*
  TObjArray binProbabilities_;      

  Double_t evaluate() const;
};

#endif
