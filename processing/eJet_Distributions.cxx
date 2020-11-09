#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPaveStats.h>
#include "Math/LorentzVector.h"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"
#include <TVector2.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TF1.h>
#include <TLatex.h>
#include <TColor.h>
#include <THStack.h>
#include <TProfile.h>
#include <iostream>
#include <fstream>
#include <TMath.h>
#include <vector>
#include <math.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <cmath>
#include <TVectorT.h>

float calc_Q_square(float inE, ROOT::Math::PtEtaPhiEVector v)
{
  return 2.*inE*v.E()*(1-TMath::Abs(TMath::Cos(v.Theta())));
}
using namespace std;
int main(int argc, char *argv[])
{
  if (argc < 2) {
    std::cout<<"Syntax: [Command] [File]"<<std::endl;
    exit(EXIT_FAILURE);
  }
  //for (int iarg = 1; iarg < argc; iarg++) {
  int iarg = 1;
  TString root_file = (TString)argv[iarg];
  
  std::cout << "Opening: " << (TString)argv[iarg] << std::endl;
  
  TFile *file = TFile::Open(root_file);
  
  if (file == NULL) {
    std::cout << " File Fail" << std::endl;
    exit(EXIT_FAILURE); 
  } 
  
  file->Print();
  
  TTree *_tree_event = dynamic_cast<TTree *>(file->Get("T"));

  if (_tree_event == NULL) {
    std::cout << " Tree Fail " << std::endl;
    exit(EXIT_FAILURE);
  }
  
  // TTreeReader Tree("T",F);
  // int nEntries = Tree.GetEntries();
  // TTreeReaderValue<int> njets(Tree,"njets");
  // TTreeReaderArray<Int_t> NConst(Tree,"nComponent");
  // TTreeReaderArray<Float_t> E(Tree,"e");
  // TTreeReaderArray<Float_t> Eta(Tree,"eta");
  // TTreeReaderArray<Float_t> Phi(Tree,"phi");
  // TTreeReaderArray<Float_t> Pt(Tree,"pt");
  
  // TTreeReaderArray<Float_t> gE(Tree,"matched_charged_truthE");
  // TTreeReaderArray<Float_t> gEta(Tree,"matched_charged_truthEta");
  // TTreeReaderArray<Float_t> gPhi(Tree,"matched_charged_truthPhi");
  // TTreeReaderArray<Float_t> gPt(Tree,"matched_charged_truthPt");
  // TTreeReaderArray<Int_t> gNComponent(Tree,"matched_charged_truthNComponent");

  // TTreeReaderValue<Float_t> electron_gE(Tree,"electron_truthE");
  // TTreeReaderValue<Float_t> electron_gEta(Tree,"electron_truthEta");
  // TTreeReaderValue<Float_t> electron_gPhi(Tree,"electron_truthPhi");
  // TTreeReaderValue<Float_t> electron_gPt(Tree,"electron_truthPt");
  
  // enum {MaxNumJets = 20,kMaxConstituents = 100};
  TString infile = argv[1];
  TFile * F = new TFile(infile);
  if (F == NULL) { std::cout << " File Fail" << std::endl; exit(EXIT_FAILURE); }
  
  TTree *T = dynamic_cast<TTree *>(F->Get("T"));
  if (_tree_event == NULL) { std::cout << " Tree Fail " << std::endl; exit(EXIT_FAILURE); }
  //TTree * T = (TTree*) F -> Get("T");
  Int_t njets;
  int nEntries = T -> GetEntries();
  const int MaxNumJets = 20;
  const int kMaxConstituents = 100;
  array<Float_t, MaxNumJets> E,Eta,Phi,Pt,gE,gEta,gPhi,gPt;
  array<Int_t, MaxNumJets> NComponent,gNComponent;
  Float_t electron_gE,electron_gEta,electron_gPhi,electron_gPt;
  array<array<Float_t, kMaxConstituents >, MaxNumJets > gComponent_Eta,gComponent_PID;

  T -> SetBranchAddress("njets",&njets);
  T -> SetBranchAddress("e",&E);
  T -> SetBranchAddress("eta",&Eta);
  T -> SetBranchAddress("phi",&Phi);
  T -> SetBranchAddress("pt",&Pt);
  T -> SetBranchAddress("nComponent",&NComponent);
  
  T -> SetBranchAddress("matched_charged_truthE",&gE);
  T -> SetBranchAddress("matched_charged_truthEta",&gEta);
  T -> SetBranchAddress("matched_charged_truthPhi",&Phi);
  T -> SetBranchAddress("matched_charged_truthPt",&Pt);
  T -> SetBranchAddress("matched_charged_truthNComponent",&gNComponent);
  T -> SetBranchAddress("matched_Constituent_truthEta", gComponent_Eta.data());
  T -> SetBranchAddress("matched_Constituent_truthPID",gComponent_PID.data());
			
  T -> SetBranchAddress("electron_truthE",&electron_gE);
  T -> SetBranchAddress("electron_truthEta",&electron_gEta);
  T -> SetBranchAddress("electron_truthPhi",&electron_gPhi);
  T -> SetBranchAddress("electron_truthPt",&electron_gPt);

  //gStyle Plotting
  gStyle->SetOptStat("emr");
  gStyle->SetStatY(0.85);
  gStyle->SetStatX(0.87);
  gStyle->SetStatW(0.15);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetHistLineWidth(2.);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes  
  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);
  
  //2D Histos
  TH2F * Tjve = new TH2F("ETrueJet_vs_Eelectron", "E^{True}_{Jet} (|#eta^{Jet}|<0.7) vs. E_{e}^{True}",100,0,25,100,0,25);
  TH2F * Rjve = new TH2F("ERecoJet_vs_Eelectron", "E^{Reco}_{Jet} (|#eta^{Jet}|<0.7) vs. E_{e}^{True}",100,0,25,100,0,25);

  //Ratio Histos
  TH1F * eoTj = new TH1F("ETrueJet_over_Eelectron", "E_{Reco}^{True}/E^{True}_{e} (|#eta^{Jet}|<0.7)",80,0,2);
  TH1F * eoRj = new TH1F("Eelectron_over_ERecoJet_over_Eelectron", "E_{Jet}^{Reco}/E^{True}_{e} (|#eta^{Jet}|<0.7)",80,0,2);
  TH1F * RjoTj = new TH1F("PRecoJet_over_ETrueJet", "P_{Jet}^{Reco} - P^{True}_{Jet} / P^{True}_{Jet}",80,-0.4,0.4);
  TH1F * justE = new TH1F("Reco Energy", "Energy",100,0,100);

  //Difference Histos
  TH1F * emTj = new TH1F("Eelectron_minus_ETrueJet", "E_{e}^{True} - E^{True}_{Jet} (|#eta^{Jet}|<0.7)",100,-20,30);
  TH1F * emRj = new TH1F("Eelectron_minus_ERecoJet", "E_{e}^{True} - E^{Reco}_{Jet} (|#eta^{Jet}|<0.7)",100,-20,30);
  //Detector Coordinate Histos
  TH1F * dPhiTj = new TH1F("dPhi_e_TrueJet", "|#Delta#varphi| (#varphi_{e} - #varphi^{True}_{Jet}) ", 32,0,M_PI);
  TH1F * dPhiRj = new TH1F("dPhi_e_RecoJet", "|#Delta#varphi| #varphi_{e} - #varphi(Jet^{Reco}_{Jet}) ", 32,0,M_PI);
  TH1F * dEtaTj = new TH1F("dEta_e_TrueJet", "|#Delta#eta| (#eta_{e} - #eta^{True}_{Jet})", 80,-10,10);
  TH1F * dEtaRj = new TH1F("dEta_e_RecoJet", "|#Delta#eta| (#eta_{e} - #eta^{Reco}_{Jet})", 80,-10,10);

  //Simple Distributions
  gStyle->SetMarkerColor(4);
  gStyle->SetLineColor(4);
  TH1F *reco_phi = new TH1F("reco_phi","Reconstructed Jet #varphi",16,-M_PI,M_PI); 
  TH1F *reco_eta = new TH1F("reco_eta","Recontsructed Jet #eta",50,-5,5);
  TH1F *reco_E = new TH1F("reco_E","Reconstructed Jet Energy",100,0,50);
  TH1F *reco_nconst = new TH1F("reco_nconst","Reconstructed Jet N Component",20,0,20);
  TH1F *reco_P = new TH1F("reco_P","Reconstructed Jet Momentum",100,0,50);
  TH1F *nconst_diff = new TH1F("nconst_diff","Truth - Reco No. Constituents",20,0,20);
  TH1F *comp_eta = new TH1F("comp_eta","Jet Component #eta",40,-4,4);
  TH1F *comp_pid = new TH1F("comp_pid","Jet Component PID",1000,-500,500);
  TH1F* Q2 = new TH1F("Q2","Q^{2}",500,0,500);
  
  gStyle->SetMarkerColor(2);
  gStyle->SetLineColor(2);
  TH1F *reco_phi_anticut = new TH1F("reco_phi(anti-cut)","Reconstructed Jet #varphi(anti-cut)",16,-M_PI,M_PI); 
  TH1F *reco_eta_anticut = new TH1F("reco_eta(anti-cut)","Recontsructed Jet #eta(anti-cut)",50,-5,5);
  TH1F *reco_E_anticut = new TH1F("reco_E(anti-cut)","Reconstructed Jet Energy(anti-cut)",100,0,50);
  TH1F *reco_nconst_anticut = new TH1F("reco_nconst(anti-cut)","Reconstructed Jet N Component(anti-cut)",20,0,20);
  TH1F *reco_P_anticut = new TH1F("reco_P(anti-cut)","Reconstructed Jet Momentum(anti-cut)",100,0,50);
  TH1F * RjoTj_anticut = new TH1F("PRecoJet_over_PTrueJet_anticut", "P_{Jet}^{True} - P_{Jet}^{True} / P^{True}_{Jet} (anti-cut) ",80,-0.4,0.4);
  TH1F *nconst_diff_anticut = new TH1F("nconst_diff_anticut","Truth - Reco No. Constituents",20,0,20);
  TH1F *comp_eta_anticut = new TH1F("comp_eta_anticut","Jet Component #eta",40,-4,4);
  TH1F *comp_pid_anticut = new TH1F("comp_pid_anticut","Jet Component PID (Anti-Cut)",1000,-500,500);
  TH1F* Q2_anticut = new TH1F("Q2_anticut","Q^{2} (anticut)",500,0,500);
  //TH1I * PID = new TH1I("PID_Histo", "PID",1000,-500,500);

  float max_DeltaR = 0.1; //reco-truth match
  int min_comp = 3;
  float minE = 4.0;
  float max_dE_E = 0.03;
  float jet_cut_counter[3] = {0};

  for (Long64_t ev = 0; ev < nEntries; ev++){
    //while ( Tree.Next() ){
    T->GetEntry(ev);
    if (ev%10000==0) fprintf(stderr,"%d: Entry %lli out of %d\n",__LINE__,ev,nEntries);
    //if (ev == 50000) break;
    //cout<<endl<<njets<<endl;
    for (int n = 0; n < njets; ++n) {//ALL truth jet loop
      //cuts
      if (std::isnan(gE[n])) continue;
      if (NComponent[n] < min_comp) continue;
      if (E[n] < minE) continue;
      if ( (abs(Eta[n]) > 1.) && (abs(Eta[n]) < 1.2) ) continue;
      // cout<<endl<<E[n]<<endl;
      // cout<<endl<<gE[n]<<endl;
      ROOT::Math::PtEtaPhiEVector Lorentz(Pt[n],Eta[n],Phi[n],E[n]);
      ROOT::Math::PtEtaPhiEVector gLorentz(gPt[n],gEta[n],gPhi[n],gE[n]);
      float dR = ROOT::Math::VectorUtil::DeltaR(Lorentz,gLorentz);
      float dE_E = (E[n] - gE[n]) / gE[n];
      if (dR > max_DeltaR) continue;

      ROOT::Math::PtEtaPhiEVector e_vector (electron_gPt,electron_gEta,electron_gPhi,electron_gE);
      Float_t Q_square = calc_Q_square(20,e_vector); //electron Beam of 20 GeV/c
      
      bool cut_to_study;
      //cut_to_study = (dR < max_DeltaR);
      cut_to_study = (abs(dE_E) < max_dE_E);
      jet_cut_counter[0]+=1.;
      if (cut_to_study)
	{
	  reco_E->Fill(E[n]);
	  reco_P->Fill(Lorentz.P());
	  reco_eta->Fill(Eta[n]);
	  reco_phi->Fill(Phi[n]);
	  reco_nconst->Fill(NComponent[n]);
	  RjoTj->Fill(dE_E);
	  jet_cut_counter[1]+=1.; //increment
	  nconst_diff->Fill(gNComponent[n] - NComponent[n]);
	  for (int i = 0; i < NComponent[n]; i++){
	    comp_eta->Fill(gComponent_Eta[n][i]);
	    comp_pid->Fill(gComponent_PID[n][i]);
	  }
	  Q2->Fill(Q_square);
	}
      else
	{
	  reco_E_anticut->Fill(E[n]);
	  reco_P_anticut->Fill(Lorentz.P());
	  reco_eta_anticut->Fill(Eta[n]);
	  reco_phi_anticut->Fill(Phi[n]);
	  reco_nconst_anticut->Fill(NComponent[n]);
	  RjoTj_anticut->Fill(dE_E);
	  nconst_diff_anticut->Fill(gNComponent[n] - NComponent[n]);
	  for (int i = 0; i < NComponent[n]; i++){
	    comp_eta_anticut->Fill(gComponent_Eta[n][i]);
	    comp_pid_anticut->Fill(gComponent_PID[n][i]);
	  }
	  Q2_anticut->Fill(Q_square);    
	}

      Float_t True_DeltaPhi = TMath::Abs(TVector2::Phi_mpi_pi(electron_gPhi - gPhi[n]));
      Float_t Reco_DeltaPhi = TMath::Abs(TVector2::Phi_mpi_pi(electron_gPhi - Phi[n]));
      dPhiTj->Fill(True_DeltaPhi);
      dPhiRj->Fill(Reco_DeltaPhi);
      dEtaTj->Fill(electron_gEta-gEta[n]);
      dEtaRj->Fill(electron_gEta-Eta[n]);
      
      //Inclusive Spectra
      Rjve->Fill(electron_gE,E[n]);
      Tjve->Fill(electron_gE,gE[n]);

      //Kinematic Cuts
      // if (True_DeltaPhi < M_PI/2) continue;
      // if (gE[n] < 3.0) continue;
      // if (TMath::Abs(eta[n]) < 0.7) continue;
      //electron/Jet Comparisons
      float eE = electron_gE;
      eoTj->Fill(gE[n]/eE);
      eoRj->Fill(E[n]/eE);
      emTj->Fill(electron_gE-gE[n]);
      emRj->Fill(electron_gE-E[n]);
    }
  }
  
  jet_cut_counter[2] = jet_cut_counter[1]/jet_cut_counter[0];
  TVectorT<double> TJet_counter(3); //total jets, jets passed, %passed
  cout<<endl<<"Fraction of Jets passing Cut:"<<jet_cut_counter[2]<<endl;
  for (int i = 0; i < 3; i++) TJet_counter[i] = jet_cut_counter[i];
  //entry loop
  
  //Write to new root file
  TFile* fout = new TFile("Histograms_Jet_Callibration.root","RECREATE");
  TJet_counter.Write("TJet_Counter");
  justE->Write();
  //PID->Write();

  //reco_E->Scale(1./reco_E->GetEntries());
  //reco_P->Scale(1./reco_P->GetEntries());
  //reco_eta->Scale(1./reco_eta->GetEntries());
  //reco_nconst->Scale(1./reco_nconst->GetEntries());
  reco_E->Write();
  reco_P->Write();
  reco_eta->Write();
  reco_phi->Write();
  reco_nconst->Write();
  RjoTj->GetXaxis()->SetTitle("dP/P");
  RjoTj->Write();
  nconst_diff->Write();
  comp_eta->Write();
  comp_pid->Write();
  
  //reco_P_anticut->Scale(1./reco_P_anticut->GetEntries());
  //reco_E_anticut->Scale(1./reco_E_anticut->GetEntries());
  //reco_eta_anticut->Scale(1./reco_eta_anticut->GetEntries());
  //reco_nconst_anticut->Scale(1./reco_nconst_anticut->GetEntries());
  reco_E_anticut->Write();
  reco_P_anticut->Write();
  reco_eta_anticut->Write();
  reco_phi_anticut->Write();
  RjoTj_anticut->Write();
  reco_nconst_anticut->Write();
  nconst_diff_anticut->Write();
  comp_eta_anticut->Write();
  comp_pid_anticut->Write();
  
  
  Tjve->GetXaxis()->SetTitle("E^{True}_{electron} [GeV]");
  Rjve->GetXaxis()->SetTitle("E^{True}_{electron} [GeV]");
  Tjve->GetYaxis()->SetTitle("E^{True}_{Jet} [GeV]");
  Rjve->GetYaxis()->SetTitle("E^{Reco}_{Jet} [GeV]");
  Tjve->Write();
  Rjve->Write();

  eoTj->GetXaxis()->SetTitle("E_{e}/E_{jet}");
  eoRj->GetXaxis()->SetTitle("E_{e}/E_{jet}");
  eoTj->Write();
  eoRj->Write();

  emTj->GetXaxis()->SetTitle("E_{e}-E_{jet}");
  emRj->GetXaxis()->SetTitle("E_{e}-E_{jet}");
  emTj->Write();
  emRj->Write();

  dPhiTj->GetXaxis()->SetTitle("#Delta#varphi");
  dPhiRj->GetXaxis()->SetTitle("#Delta#varphi");
  dPhiTj->Write();
  dPhiRj->Write();

  dEtaTj->GetXaxis()->SetTitle("#Delta#eta");
  dEtaRj->GetXaxis()->SetTitle("#Delta#eta");
  dEtaTj->Write();
  dEtaRj->Write();

  Q2->Write();
  Q2_anticut->Write();
  // H_dR->Write();
  // H_NExtra_Matches->Write();

  fout->Close();
  
}

