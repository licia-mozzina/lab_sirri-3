#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "TCanvas.h"
#include "TH1D.h"
#include "TRandom.h"

#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"
//#include "RooUnfoldSvd.h"
//#include "RooUnfoldTUnfold.h"
//#include "RooUnfoldIds.h"
#endif

// Dummy variable for rejected data
const Double_t cutdummy = -99999.0;

// Gaussian smearing, systematic translation, and variable inefficiency
Double_t smear(Double_t xt) {
  Double_t xeff = 0.3 + (1.0 - 0.3) / 20 * (xt + 10.0); // efficiency
  Double_t x = gRandom->Rndm();
  if (x > xeff)
    return cutdummy;
  Double_t xsmear = gRandom->Gaus(-2.5, 0.2); // bias and smear
  return xt + xsmear;
}

void unfold_psexp() {
  // gSystem->Load("~/RooUnfold/build/libRooUnfold.so");

  cout << "==================================== TRAIN "
          "===================================="
       << endl;
  RooUnfoldResponse response(40, -10.0, 10.0);
  RooUnfoldResponse response_psexp(40, -10.0, 10.0);   // For the pseudo experiments part

  // Train with a Breit-Wigner, mean 0.3 and width 2.5.
  for (Int_t i = 0; i < 100000; i++) {
    Double_t xt = gRandom->BreitWigner(0.3, 2.5);
    Double_t x = smear(xt);
    if (x != cutdummy) {
      response.Fill(x, xt);
      response_psexp.Fill(x, xt);
    } else {
      response.Miss(xt);
      response_psexp.Miss(xt);
    }
  }

  cout << "==================================== TEST "
          "====================================="
       << endl;
  TH1D *hTrue = new TH1D("true", "Test Truth", 40, -10.0, 10.0);
  TH1D *hMeas = new TH1D("meas", "Test Measured", 40, -10.0, 10.0);

  // Test with a Gaussian, mean 0 and width 2.
  for (Int_t i = 0; i < 10000; i++) {
    Double_t xt = gRandom->Gaus(0.0, 2.0), x = smear(xt);
    hTrue->Fill(xt);
    if (x != cutdummy)
      hMeas->Fill(x);
  }

  cout << "==================================== UNFOLD "
          "==================================="
       << endl;
  RooUnfoldBayes unfold(&response, hMeas, 4); // OR
  // RooUnfoldSvd     unfold (&response, hMeas, 20);   // OR
  // RooUnfoldTUnfold unfold (&response, hMeas);       // OR
  // RooUnfoldIds     unfold (&response, hMeas, 1);

  TH1D *hUnfold = (TH1D *)unfold.Hunfold();

  //****************************************************************

  // pseudoexperiments part

  // Printing statistical uncertainties estimated via the D'Agostini's method
  std::vector<Double_t> cov_ij_dag;

  for (Int_t i = 1; i < hUnfold->GetNbinsX(); i++) {
    std::cout << "Uncertainties via D'Agostini's method" << hUnfold->GetBinError(i) << endl;
    cov_ij_dag.push_back(hUnfold->GetBinError(i));
  }

  // Histograms for pseudo experiments
  TH1D *hTrue_psexp = new TH1D("h_truth", "Test Truth", 40, -10.0, 10.0);
  TH1D *hMeas_psexp = new TH1D("h_meas", "Test Measured", 40, -10.0, 10.0);


  Int_t nToy = 100;  //number of pseudo experiments

  // Array of histograms for each unfolded pseudo experiment
  TH1D *hUnfold_psexp[nToy];
  char name[20];
  char title[100];
  for (Int_t i = 0; i < nToy; i++) {
    sprintf(name, "h_psexp_%d", i);
    sprintf(title, "title of h_pesexp_%d", i);
    hUnfold_psexp[i] = new TH1D(name, title, 40, -10.0, 10.0);
  }


  std::vector<Double_t> bin_content_mean;

  // Pseudo experiments loop
  for (Int_t k = 0; k < nToy; k++) {

    // Preparing the "measured" spectrum such that event counts in each bin are drawn from Poisson distribution
    for (Int_t i = 0; i < hMeas->GetNbinsX(); i++) {
      Double_t xt = gRandom->Poisson(hMeas->GetBinContent(i)),
               x = smear(xt); // devo settarlo come valore al contenuto del bin
      hTrue_psexp->Fill(xt);
      if (x != cutdummy)
        hMeas_psexp->Fill(x);
    }

    // Unfolding
    RooUnfoldBayes unfold_psexp(&response_psexp, hMeas_psexp, 4);
    hUnfold_psexp[k] = (TH1D *)unfold_psexp.Hunfold();

    for (Int_t i = 1; i < hUnfold_psexp[k]->GetNbinsX() + 1; ++i) {
      if (k == 0) {
        bin_content_mean.push_back(hUnfold_psexp[k]->GetBinContent(i) / nToy); // from first pesudo exp we set the number of vector elements
                                                                       // we are already averaging over nToy number of pseudo experiments     
      } else {
        bin_content_mean[i] += (hUnfold_psexp[k]->GetBinContent(i) / nToy); // sum bin values for each pseudo exp
      }
    }
  }


  // Estimating the covariance matrix for the pseudo experiments
  std::vector<Double_t> cov_ij_psexp;
  std::vector<Double_t> cov_ii_psexp;

  // Looping to compute the uncertainties
  for (Int_t k = 0; k < nToy; ++k) {
    if (k == 0) { // from first pesudo exp we set the number of vector elements
      for (Int_t ibin = 1; ibin < hUnfold_psexp[k]->GetNbinsX() + 1; ibin++) {
        for (Int_t jbin = 1; jbin < hUnfold_psexp[k]->GetNbinsX() + 1; jbin++) {

          Double_t cov_ij_N =
              (hUnfold_psexp[k]->GetBinContent(ibin) - bin_content_mean[ibin]) *
              (hUnfold_psexp[k]->GetBinContent(jbin) - bin_content_mean[jbin]);

          cov_ij_psexp.push_back(cov_ij_N / nToy); // averaging over all pseudo experiments
          
          if (ibin == jbin) cov_ii_psexp.push_back(cov_ij_N / nToy);
        }
      }
    } else {
      for (Int_t ibin = 1; ibin < hUnfold_psexp[k]->GetNbinsX() + 1; ibin++) {
        for (Int_t jbin = 1; jbin < hUnfold_psexp[k]->GetNbinsX() + 1; jbin++) {
          Int_t index = ibin * hUnfold_psexp[k]->GetNbinsX() + jbin;

          Double_t cov_ij_N =
              (hUnfold_psexp[k]->GetBinContent(ibin) - bin_content_mean[ibin]) *
              (hUnfold_psexp[k]->GetBinContent(jbin) - bin_content_mean[jbin]);

          cov_ij_psexp[index] += cov_ij_N / nToy;

          if (ibin == jbin) cov_ii_psexp[index] += cov_ij_N / nToy;
        }
      }
    }
  }

  // Comparing the uncertainties
  for (Int_t i = 0; i < cov_ii_psexp.size(); ++i) {
    Double_t cov_ij_psexp_root = TMath::Sqrt(cov_ii_psexp[i]);
    printf("D'Agostini %f | Pseudo experiment %f \n", cov_ij_dag[i],
           cov_ij_psexp_root);
  } // fanno un po' schifo
}

#ifndef __CINT__
int main() {
  unfold_psexp();
  return 0;
} // Main program when run stand-alone
#endif
