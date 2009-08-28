void bgEstTemplateShapeBias()
{
  TString dqmDirectoryName = "DQMData/harvested/TTplusJets/zMuTauAnalyzer/";

  TObjArray dqmSubDirectoryNames;
  dqmSubDirectoryNames.Add(new TObjString("afterEvtSelDiTauCandidateForMuTauAcoplanarity12_beforeEvtSelDiTauCandidateForMuTauMt1MET/"));
  dqmSubDirectoryNames.Add(new TObjString("afterEvtSelDiTauCandidateForMuTauMt1MET_beforeEvtSelDiTauCandidateForMuTauPzetaDiff/"));
  dqmSubDirectoryNames.Add(new TObjString("afterEvtSelDiTauCandidateForMuTauPzetaDiff/"));

  TObjArray labels;
  labels.Add(new TObjString("before M_{T}^{#mu + MET} Cut"));
  labels.Add(new TObjString("after M_{T}^{#mu + MET} Cut, before Cut on P_{#zeta} - 1.5*P_{#zeta}^{vis}"));
  labels.Add(new TObjString("after Cut on P_{#zeta} - 1.5*P_{#zeta}^{vis}"));

  TObjArray dqmMonitorElementNames;
  dqmMonitorElementNames.Add(new TObjString("MuonQuantities/MuonPt"));
  dqmMonitorElementNames.Add(new TObjString("MuonQuantities/MuonEta"));
  dqmMonitorElementNames.Add(new TObjString("TauQuantities/TauPt"));
  dqmMonitorElementNames.Add(new TObjString("TauQuantities/TauEta"));
  dqmMonitorElementNames.Add(new TObjString("DiTauCandidateQuantities/DPhi12"));
  dqmMonitorElementNames.Add(new TObjString("DiTauCandidateQuantities/DR12"));
  dqmMonitorElementNames.Add(new TObjString("DiTauCandidateQuantities/DPhi1MET"));
  dqmMonitorElementNames.Add(new TObjString("DiTauCandidateQuantities/DPhi2MET"));
  dqmMonitorElementNames.Add(new TObjString("DiTauCandidateQuantities/PzetaDiff"));
  dqmMonitorElementNames.Add(new TObjString("DiTauCandidateQuantities/VisMass"));

  TString outputFileName_unnormalized = "bgEstTemplateShapeBias_unnormalized.ps";
  TString outputFileName_normalized = "bgEstTemplateShapeBias_normalized.ps";

  TString inputFileName = "../../Configuration/test/plotsZtoMuTau_all_shrinkingCone.root";

  showTemplateShapeBiasAll(inputFileName, dqmDirectoryName, dqmSubDirectoryNames, labels, dqmMonitorElementNames, 
			   false, outputFileName_unnormalized);
  showTemplateShapeBiasAll(inputFileName, dqmDirectoryName, dqmSubDirectoryNames, labels, dqmMonitorElementNames, 
			   true, outputFileName_normalized);
}

void showTemplateShapeBiasAll(const TString& inputFileName, const TString& dqmDirectoryName, 
			      const TObjArray& dqmSubDirectoryNames, const TObjArray& labels, const TObjArray& dqmMonitorElementNames,
			      bool normalize, const TString& outputFileName)
{
  if ( dqmSubDirectoryNames.GetEntries() != labels.GetEntries() ) {
    std::cerr << "Error in <showTemplateShapeBiasAll>: size of dqmSubDirectoryNames and labels collections don't match" << 
	      << " --> histograms will NOT be plotted !!" << std::endl;
    return;
  }

  TCanvas* canvas = new TCanvas("canvas", "canvas", 1, 1, 800, 600);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);

  TPostScript* ps = new TPostScript(outputFileName, 112);
  
  TFile inputFile(inputFileName);

  int numPlots = dqmMonitorElementNames.GetEntries();
  for ( int iPlot = 0; iPlot < numPlots; ++iPlot ) {
    TString dqmMonitorElementName = ((TObjString*)dqmMonitorElementNames.At(iPlot))->GetString();

    TObjArray dqmMonitorElements;
    
    int numDistributions = dqmSubDirectoryNames.GetEntries();
    for ( int iDistribution = 0; iDistribution < numDistributions; ++iDistribution ) {
      TString dqmSubDirectoryName = ((TObjString*)dqmSubDirectoryNames.At(iDistribution))->GetString();

      TString dqmMonitorElementName_full = TString(dqmDirectoryName).Append(dqmSubDirectoryName).Append(dqmMonitorElementName);
      std::cout << "dqmMonitorElementName_full = " << dqmMonitorElementName_full << std::endl;

      TH1* dqmMonitorElement = (TH1*)inputFile.Get(dqmMonitorElementName_full);
      std::cout << "dqmMonitorElement = " << dqmMonitorElement <<  std::endl;

      TH1* dqmMonitorElement_cloned = (TH1*)dqmMonitorElement->Clone();
      
      if ( normalize ) dqmMonitorElement_cloned->Scale(1./dqmMonitorElement_cloned->Integral());
      
      dqmMonitorElements.Add(dqmMonitorElement_cloned);
    }

    showTemplateShapeBiasOne(canvas, ps, dqmMonitorElements, labels);
  }

  delete ps;

  delete canvas;
}

void showTemplateShapeBiasOne(TCanvas* canvas, TPostScript* ps, const TObjArray& dqmMonitorElements, const TObjArray& labels)
{
  if ( dqmMonitorElements.GetEntries() != labels.GetEntries() ) {
    std::cerr << "Error in <showTemplateShapeBiasOne>: size of dqmMonitorElements and labels collections don't match" << 
	      << " --> histograms will NOT be plotted !!" << std::endl;
    return;
  }

//--- determine scale of y-axis
  Float_t yMax = 0.;

  int numDistributions = dqmMonitorElements.GetEntries();
  for ( int iDistribution = 0; iDistribution < numDistributions; ++iDistribution ) {
    TH1* dqmMonitorElement = (TH1*)dqmMonitorElements.At(iDistribution);

    if ( dqmMonitorElement->GetMaximum() > yMax ) yMax = dqmMonitorElement->GetMaximum();
  }

  TLegend legend(0.53, 0.73, 0.88, 0.89);
  legend.SetBorderSize(0);
  legend.SetFillColor(0);

  for ( int iDistribution = 0; iDistribution < numDistributions; ++iDistribution ) {
    TH1* dqmMonitorElement = (TH1*)dqmMonitorElements.At(iDistribution);

    dqmMonitorElement->SetMaximum(1.3*yMax);
    dqmMonitorElement->SetMinimum(0.);

    dqmMonitorElement->SetStats(false);

    dqmMonitorElement->SetLineColor(2 + iDistribution);
    dqmMonitorElement->SetLineWidth(2);

    const char* drawOption = ( iDistribution == 0 ) ? "" : "same";
    dqmMonitorElement->Draw(drawOption);

    TString label = ((TObjString*)labels.At(iDistribution))->GetString();
    std::cout << "label = " << label << std::endl;

    legend.AddEntry(dqmMonitorElement, label, "l");
  }

  legend.Draw();

  canvas->Update();
  
  ps->NewPage();
}
