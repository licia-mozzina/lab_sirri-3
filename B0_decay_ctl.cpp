#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooFitResult.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"
#include "TH2.h"

using namespace RooFit;



int B0_decay_ctl() {
  // first declare the observable
  RooRealVar x("x", "m(p-pbar)", 5090., 5590.);
  x.setBins(24);

  // background model
  RooRealVar c("c", "exponential par", -0.001, -0.01, 0.); // mancavano i punti decimali, era confuso
  RooExponential bkg("bkg", "Comb. bkg.", x, c);

  // first gaussian
  RooRealVar m0("m0", "B0 mass", 5279., 5220., 5320.);
  RooRealVar s0("s0", "B0 width", 10., 2., 50.);
  RooGaussian gaus0("gaus0", "B0 peak", x, m0, s0);

  // second gaussian
  RooRealVar m1("m1", "B0s mass", 5380., 5320., 5420.);
  RooRealVar s1("s1", "B0s width", 10., 2., 50.);
  RooGaussian gaus1("gaus1", "B0s peak", x, m1, s1);

  // Background in control region
  RooRealVar bg("bg", "bg_control_region", 4000., 5000.);
  RooExponential model_ctl("model_ctl", "control_bkg", x, c);
  c.setVal(-1.0e-3);

  /*
  // composite model
  // not extended model  3 p.d.s. => 2 fractions
  //   f0 * gaus0(x) + f1 * gaus1(x) + (1 -f0 -f1)* bkg(x)

  RooRealVar f0("f0", "fraction of B0" , 0.2, 0.00001, 1);
  RooRealVar f1("f1", "fraction of B0s", 0.2, 0.00001, 1);

  RooAddPdf model("model", "model",
     RooArgList(gaus0, gaus1, bkg), RooArgList(f0, f1) );
  */

  // composite model
  // extended model  3 p.d.s. => 3 expected n. of events
  //   N0 * gaus0(x) + N1 * gaus1(x) + Nb* bkg(x)

  RooRealVar n0("n0", "number of B0 events", 100, 0, 10000);
  RooRealVar n1("n1", "number of B0s events", 10, 0, 10000);
  RooRealVar nb("nb", "number of background events", 1000, 0, 100000);

  RooAddPdf model("extended model", "extended model", RooArgList(gaus0, gaus1, bkg), RooArgList(n0, n1, nb));

  // Reading data from input file
  RooDataSet *data = RooDataSet::read("rarest_b0_decay.dat", x, "x");

  // Generating data for control region
  RooDataSet *data_ctl = model_ctl.generate(x, 1e4); 

  // Define category to distinguish physics and control region 
  // and constructing a joint data sample and pdf
  RooCategory sample("sample", "Sample category", {
                  {"Physics", 0},
                  {"Control", 1}
                });

  RooDataSet combData("combData", "combined data", RooArgSet(x, bg), Index(sample), Import({{"Physics", data}, {"Control", data_ctl}}));
  RooSimultaneous simPdf("simPdf", "simultaneous pdf", {{"Physics", &model}, {"Control", &model_ctl}}, sample);

  // Simultaneous fit
  RooFitResult* comb_fit = simPdf.fitTo(combData, Save());
  comb_fit->Print();

  // Plotting physics sample
  RooPlot *frame1 = x.frame(Title("Physics sample"));
  combData.plotOn(frame1, Cut("sample==sample::Physics"));
  frame1->GetXaxis()->SetTitle("m(p#bar{p})[MeV/c^{2}]");

  // Plotting "physics" slice of simultaneous pdf
  simPdf.plotOn(frame1, Slice(sample, "Physics"), ProjWData(sample, combData));
  simPdf.plotOn(frame1, Slice(sample, "Physics"), Components("bkg"), ProjWData(sample, combData), LineStyle(kDashed));

  // Plotting "sample" slice of simultaneous pdf
  RooPlot *frame2 = x.frame(Title("Control sample")); 
  RooAbsData* slicedData{combData.reduce(Cut("sample==sample::Control"))};
  slicedData->plotOn(frame2);
  simPdf.plotOn(frame2, ProjWData(sample, *slicedData));
  simPdf.plotOn(frame2, Components("model_ctl"), ProjWData(sample, *slicedData), LineStyle(kDashed));
  frame2->GetXaxis()->SetTitle("m(p#bar{p})[MeV/c^{2}]");
 
  // Plotting both samples
  RooPlot *frame3 = x.frame(Bins(24), Title("Both samples"));
  combData.plotOn(frame3);
  simPdf.plotOn(frame3, ProjWData(sample, combData));
  simPdf.plotOn(frame3, Components("bkg, model_ctl"), ProjWData(sample, combData),
                LineStyle(kDashed));
  frame3->GetXaxis()->SetTitle("m(p#bar{p})[MeV/c^{2}]");

  // Creating canvas and saving plot 
  TCanvas *c1 = new TCanvas("B0_decay_ctl", "B0 decay control model", 1200, 400);
  c1->Divide(3);
  auto draw = [&](int i, RooPlot & frame) {
      c1->cd(i);
      gPad->SetLeftMargin(0.15);
      frame.GetYaxis()->SetTitleOffset(1.4);
      frame.Draw();
  };
  draw(1, *frame1);
  draw(2, *frame2);
  draw(3, *frame3);

  c1->Print("B0_decay_ctl.png");
  
  return 0;
}