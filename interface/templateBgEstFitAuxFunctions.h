#ifndef TauAnalysis_BgEstimationTools_templateBgEstFitAuxFunctions_h
#define TauAnalysis_BgEstimationTools_templateBgEstFitAuxFunctions_h

#include <TH1.h>
#include <TArrayD.h>
#include <TVectorD.h>
#include <TMatrixD.h>

#include <RooFitResult.h>

enum { kCoherent, kIncoherent };

double getSampledPull(double, double, double);

void sampleHistogram_stat(const TH1*, TH1*);
void sampleHistogram_sys(TH1*, const TH1*, double, double, double, int);

TArrayD getBinning(const TH1*);

double getIntegral(const TH1*);
double getIntegral(const TH1*, double, double);

void makeHistogramPositive(TH1*);

void unpackFitResult(const RooFitResult*, TVectorD&, TMatrixD&);

#endif
