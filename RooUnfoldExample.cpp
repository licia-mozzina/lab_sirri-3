//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id$
//
// Description:
//      Simple example usage of the RooUnfold package using toy MC.
//
// Authors: Tim Adye <T.J.Adye@rl.ac.uk> and Fergus Wilson <fwilson@slac.stanford.edu>
//
//==============================================================================

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"
#include "TCanvas.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
//#include "RooUnfoldSvd.h"
//#include "RooUnfoldTUnfold.h"
//#include "RooUnfoldIds.h"
#endif


//==============================================================================
// Global definitions
//==============================================================================

const Double_t cutdummy= -99999.0;

//==============================================================================
// Gaussian smearing, systematic translation, and variable inefficiency
//==============================================================================

Double_t smear (Double_t xt)
{
  Double_t xeff= 0.3 + (1.0-0.3)/20*(xt+10.0);  // efficiency
  Double_t x= gRandom->Rndm();
  if (x>xeff) return cutdummy;
  Double_t xsmear= gRandom->Gaus(-2.5,0.2);     // bias and smear
  return xt+xsmear;
}

//==============================================================================
// Example Unfolding
//==============================================================================

void RooUnfoldExample()
{  
  //gSystem->Load("~/RooUnfold/build/libRooUnfold.so");

  cout << "==================================== TRAIN ====================================" << endl;
  RooUnfoldResponse response (40, -10.0, 10.0);

  // Train with a Breit-Wigner, mean 0.3 and width 2.5.
  for (Int_t i= 0; i<100000; i++) {
    Double_t xt= gRandom->BreitWigner (0.3, 2.5);
    Double_t x= smear (xt);
    if (x!=cutdummy)
      response.Fill (x, xt);
    else
      response.Miss (xt);
  }

  cout << "==================================== TEST =====================================" << endl;
  TH1D* hTrue= new TH1D ("true", "Test Truth",    40, -10.0, 10.0);
  TH1D* hMeas= new TH1D ("meas", "Test Measured", 40, -10.0, 10.0);
  // Test with a Gaussian, mean 0 and width 2.
  for (Int_t i=0; i<10000; i++) {
    Double_t xt= gRandom->Gaus (0.0, 2.0), x= smear (xt);
    hTrue->Fill(xt);
    if (x!=cutdummy) hMeas->Fill(x);
  }

  cout << "==================================== UNFOLD ===================================" << endl;
  RooUnfoldBayes   unfold (&response, hMeas, 4);    // OR
//RooUnfoldSvd     unfold (&response, hMeas, 20);   // OR
//RooUnfoldTUnfold unfold (&response, hMeas);       // OR
//RooUnfoldIds     unfold (&response, hMeas, 1);

  TH1D* hUnfold= (TH1D*) unfold.Hunfold();

  TCanvas* c1= new TCanvas("c1","Truth vs measured", 1600, 800);

  unfold.PrintTable (cout, hTrue);
  hUnfold->Draw();
  hMeas->Draw("SAME");
  hTrue->SetLineColor(8);
  hTrue->Draw("SAME");

  c1->SaveAs("RooUnfoldExample.pdf");


//****************************************************************
  
  TMatrixD responseMatrix = response.Mresponse();
  
  TCanvas* c2 = new TCanvas("c2","response matrix", 1600, 800);
  responseMatrix.Draw("colz");

  auto responseInv = *static_cast<TMatrixD*>(responseMatrix.Clone());
  responseInv.Invert();

  TCanvas* c3 = new TCanvas("c3","inverse response matrix", 1600, 800);
  responseInv.Draw("colz");

  auto model_meas_Rinv = static_cast<TH1D*>(hMeas->Clone("model_meas_Rinv"));
  model_meas_Rinv->Reset();

  for (Int_t ibin = 1; ibin < model_meas_Rinv->GetNbinsX() + 1; ibin++) {
    Double_t matsum = 0.;
    for (Int_t jbin = 1; jbin < model_meas_Rinv->GetNbinsX() + 1; jbin++) {
      matsum += hMeas->GetBinContent(jbin) * responseInv[ibin - 1][jbin - 1];
    }
    model_meas_Rinv->SetBinContent(ibin, matsum);
  }

  TCanvas* c4 = new TCanvas("c4","matrix inversion approach", 1600, 800);
  c4->Divide(2,2);
  c4->cd(1);
  hTrue->Draw("HIST");
  c4->cd(2);
  hMeas->Draw("PE SAME");
  c4->cd(3);
  model_meas_Rinv->SetLineColor(kViolet);
  model_meas_Rinv->Draw("HIST SAME");
  c4->cd(4);
  model_meas_Rinv->Draw("HIST");

//*****************************************************************


//pseudoexperiments part

for (Int_t i = 1; i < hUnfold.GetNbinsX() + 1; i++) { //perchè partiamo da uno?
  print hUnfold->GetBinError(i);
}

TH1D* hTrue_psexp= new TH1D ("h_truth", "Test Truth",  40, -10.0, 10.0);
TH1D* hMeas_psexp= new TH1D ("h_meas", "Test Measured", 40, -10.0, 10.0);

// pseudo epx loop
Int_t nToy = 100;
TTree * t = new TTree("t", "pseudoexp_spectra");
RooUnfoldResponse response_psexp (40, -10.0, 10.0);
RooUnfoldBayes pseudoexp_spectrum (&response_psexp, hMeas, 4);
t->Branch("pseudoexp_spectrum", &pseudoexp_spectrum);

for (Int_t j = 0; j < nToy; j++) {
  
  // RooUnfoldResponse response_psexp (40, -10.0, 10.0);
  // RooUnfoldBayes   unfold_psexp (&response_psexp, hMeas, 4);    
  TH1D* hUnfold_psexp = (TH1D*) unfold_psexp.Hunfold();
  
  for (Int_t i=0; i< hMeas->GetNBins(); i++) {
    Double_t xt= gRandom->Poisson(hMeas->GetBinContent(i)), x= smear (xt); //devo settarlo come valore al contenuto del bin
    hTrue_psexp->Fill(xt);
    t->Fill();
    if (x!=cutdummy) hMeas_psexp->Fill(x);
  }
}





}

#ifndef __CINT__
int main () { RooUnfoldExample(); return 0; }  // Main program when run stand-alone
#endif
