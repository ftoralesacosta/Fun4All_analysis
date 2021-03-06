#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPaveStats.h>
#include "Math/LorentzVector.h"
#include "Math/Vector4D.h"
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

float calc_Q_square(float inE, TLorentzVector v)
{
  return 2.*inE*v.E()*(1-TMath::Abs(v.CosTheta()));
}

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

  enum {MaxNumJets = 20,kMaxConstituents = 100};

  //Declare Leaf Types
  Int_t           event;
  Int_t           njets;
  Int_t           nAlltruthjets;
  
  std::array <Int_t,MaxNumJets>           id;
  std::array <Int_t,MaxNumJets>           nComponent;
  std::array <Float_t,MaxNumJets>         eta;
  std::array <Float_t,MaxNumJets>         phi;
  std::array <Float_t,MaxNumJets>         e;
  std::array <Float_t,MaxNumJets>         pt;

  std::array <Int_t,MaxNumJets>           matched_truthID;
  std::array <Int_t,MaxNumJets>           matched_truthNComponent;
  std::array <Float_t,MaxNumJets>         matched_truthEta;
  std::array <Float_t,MaxNumJets>         matched_truthPhi;
  std::array <Float_t,MaxNumJets>         matched_truthE;
  std::array <Float_t,MaxNumJets>         matched_truthPt;

  std::array<std::array<Float_t, kMaxConstituents >, MaxNumJets > matched_Constituent_truthPID;
  std::array<std::array<Float_t, kMaxConstituents >, MaxNumJets > matched_Constituent_truthCharge;
  std::array<std::array<Float_t, kMaxConstituents >, MaxNumJets > matched_Constituent_truthEta;
  std::array<std::array<Float_t, kMaxConstituents >, MaxNumJets > matched_Constituent_truthE;
  std::array<std::array<Float_t, kMaxConstituents >, MaxNumJets > matched_Constituent_truthPhi;
  std::array<std::array<Float_t, kMaxConstituents >, MaxNumJets > matched_Constituent_truthPt;

  std::array <Int_t,MaxNumJets>           all_truthID;
  std::array <Int_t,MaxNumJets>           all_truthNComponent;
  std::array <Float_t,MaxNumJets>         all_truthEta;
  std::array <Float_t,MaxNumJets>         all_truthPhi;
  std::array <Float_t,MaxNumJets>         all_truthE;
  std::array <Float_t,MaxNumJets>         all_truthPt;

  std::array<std::array<Float_t, kMaxConstituents >, MaxNumJets > all_Constituent_truthPID;
  std::array<std::array<Float_t, kMaxConstituents >, MaxNumJets > all_Constituent_truthCharge;
  std::array<std::array<Float_t, kMaxConstituents >, MaxNumJets > all_Constituent_truthEta;
  std::array<std::array<Float_t, kMaxConstituents >, MaxNumJets > all_Constituent_truthE;
  std::array<std::array<Float_t, kMaxConstituents >, MaxNumJets > all_Constituent_truthPhi;
  std::array<std::array<Float_t, kMaxConstituents >, MaxNumJets > all_Constituent_truthPt;
  

  Float_t         electron_truthEta;
  Float_t         electron_truthPhi;
  Float_t         electron_truthE;
  Float_t         electron_truthPt;
  Float_t         electron_truthpX;
  Float_t         electron_truthpY;
  Float_t         electron_truthpZ;
  Int_t           electron_truthPID;
  Int_t           electron_truthParentID;

  //Declare Branches

  TBranch        *b_event;
  TBranch        *b_njets;
  TBranch        *b_nAlltruthjets;

  //Reco
  TBranch        *b_id;
  TBranch        *b_nComponent;
  TBranch        *b_eta;
  TBranch        *b_phi;
  TBranch        *b_e;
  TBranch        *b_pt;

  //Truth
  TBranch        *b_matched_truthID;
  TBranch        *b_matched_truthNComponent;
  TBranch        *b_matched_truthEta;
  TBranch        *b_matched_truthPhi;
  TBranch        *b_matched_truthE;
  TBranch        *b_matched_truthPt;

  TBranch        *b_matched_Constituent_truthPID;
  TBranch        *b_matched_Constituent_truthCharge;
  TBranch        *b_matched_Constituent_truthNComponent;
  TBranch        *b_matched_Constituent_truthEta;
  TBranch        *b_matched_Constituent_truthPhi;
  TBranch        *b_matched_Constituent_truthE;
  TBranch        *b_matched_Constituent_truthPt;

  TBranch        *b_all_truthID;
  TBranch        *b_all_truthNComponent;
  TBranch        *b_all_truthEta;
  TBranch        *b_all_truthPhi;
  TBranch        *b_all_truthE;
  TBranch        *b_all_truthPt;

  TBranch        *b_all_Constituent_truthPID;
  TBranch        *b_all_Constituent_truthCharge;
  TBranch        *b_all_Constituent_truthNComponent;
  TBranch        *b_all_Constituent_truthEta;
  TBranch        *b_all_Constituent_truthPhi;
  TBranch        *b_all_Constituent_truthE;
  TBranch        *b_all_Constituent_truthPt;  
  
  TBranch        *b_electron_truthEta;
  TBranch        *b_electron_truthPhi;
  TBranch        *b_electron_truthE;
  TBranch        *b_electron_truthPt;
  TBranch        *b_electron_truthpX;
  TBranch        *b_electron_truthpY;
  TBranch        *b_electron_truthpZ;
  TBranch        *b_electron_truthPID;
  TBranch        *b_electron_truthParentID;

  //Set Branch Addresses
  _tree_event->SetBranchAddress("event", &event, &b_event);
  _tree_event->SetBranchAddress("njets",&njets, &b_njets);
  _tree_event->SetBranchAddress("nAlltruthjets",&nAlltruthjets, &b_nAlltruthjets);

  _tree_event->SetBranchAddress("id", &id, &b_id);
  _tree_event->SetBranchAddress("nComponent", &nComponent, &b_nComponent);
  _tree_event->SetBranchAddress("eta", &eta, &b_eta);
  _tree_event->SetBranchAddress("phi", &phi, &b_phi);
  _tree_event->SetBranchAddress("e", &e, &b_e);
  _tree_event->SetBranchAddress("pt", &pt, &b_pt);


  _tree_event->SetBranchAddress("matched_truthID", &matched_truthID, &b_matched_truthID);
  _tree_event->SetBranchAddress("matched_truthNComponent", &matched_truthNComponent, &b_matched_truthNComponent);
  _tree_event->SetBranchAddress("matched_truthEta", &matched_truthEta, &b_matched_truthEta);
  _tree_event->SetBranchAddress("matched_truthPhi", &matched_truthPhi, &b_matched_truthPhi);
  _tree_event->SetBranchAddress("matched_truthE", &matched_truthE, &b_matched_truthE);
  _tree_event->SetBranchAddress("matched_truthPt", &matched_truthPt, &b_matched_truthPt);

  _tree_event->SetBranchAddress("matched_Constituent_truthPID", matched_Constituent_truthPID.data(), &b_matched_Constituent_truthPID);
  _tree_event->SetBranchAddress("matched_Constituent_truthCharge", matched_Constituent_truthCharge.data(), &b_matched_Constituent_truthCharge);
  _tree_event->SetBranchAddress("matched_Constituent_truthEta", matched_Constituent_truthEta.data(), &b_matched_Constituent_truthEta);
  _tree_event->SetBranchAddress("matched_Constituent_truthE", matched_Constituent_truthE.data(), &b_matched_Constituent_truthE);
  _tree_event->SetBranchAddress("matched_Constituent_truthPhi", matched_Constituent_truthPhi.data(), &b_matched_Constituent_truthPhi);
  _tree_event->SetBranchAddress("matched_Constituent_truthPt", matched_Constituent_truthPt.data(), &b_matched_Constituent_truthPt);


  _tree_event->SetBranchAddress("all_truthID", &all_truthID, &b_all_truthID);
  _tree_event->SetBranchAddress("all_truthNComponent", &all_truthNComponent, &b_all_truthNComponent);
  _tree_event->SetBranchAddress("all_truthEta", &all_truthEta, &b_all_truthEta);
  _tree_event->SetBranchAddress("all_truthPhi", &all_truthPhi, &b_all_truthPhi);
  _tree_event->SetBranchAddress("all_truthE", &all_truthE, &b_all_truthE);
  _tree_event->SetBranchAddress("all_truthPt", &all_truthPt, &b_all_truthPt);
 
  _tree_event->SetBranchAddress("all_Constituent_truthPID", all_Constituent_truthPID.data(), &b_all_Constituent_truthPID);
  _tree_event->SetBranchAddress("all_Constituent_truthCharge", all_Constituent_truthCharge.data(), &b_all_Constituent_truthCharge);
  _tree_event->SetBranchAddress("all_Constituent_truthEta", all_Constituent_truthEta.data(), &b_all_Constituent_truthEta);
  _tree_event->SetBranchAddress("all_Constituent_truthE", all_Constituent_truthE.data(), &b_all_Constituent_truthE);
  _tree_event->SetBranchAddress("all_Constituent_truthPhi", all_Constituent_truthPhi.data(), &b_all_Constituent_truthPhi);
  _tree_event->SetBranchAddress("all_Constituent_truthPt", all_Constituent_truthPt.data(), &b_all_Constituent_truthPt);

  
  _tree_event->SetBranchAddress("electron_truthEta", &electron_truthEta, &b_electron_truthEta);
  _tree_event->SetBranchAddress("electron_truthPhi", &electron_truthPhi, &b_electron_truthPhi);
  _tree_event->SetBranchAddress("electron_truthE", &electron_truthE, &b_electron_truthE);
  _tree_event->SetBranchAddress("electron_truthPt", &electron_truthPt, &b_electron_truthPt);
  _tree_event->SetBranchAddress("electron_truthpX", &electron_truthpX, &b_electron_truthpX);
  _tree_event->SetBranchAddress("electron_truthpY", &electron_truthpY, &b_electron_truthpY);
  _tree_event->SetBranchAddress("electron_truthpZ", &electron_truthpZ, &b_electron_truthpZ);
  _tree_event->SetBranchAddress("electron_truthPt", &electron_truthPt, &b_electron_truthPt);
  _tree_event->SetBranchAddress("electron_truthPID", &electron_truthPID, &b_electron_truthPID);

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
  TH1F * RjoTj = new TH1F("ERecoJet_over_ETrueJet", "E_{Jet}^{True}/E^{Reco}_{Jet} (|#eta^{Jet}|<0.7)",80,0,2);
  TH1F * justE = new TH1F("Reco Energy", "Energy",100,0,100);

  //Difference Histos
  TH1F * emTj = new TH1F("Eelectron_minus_ETrueJet", "E_{e}^{True} - E^{True}_{Jet} (|#eta^{Jet}|<0.7)",100,-20,30);
  TH1F * emRj = new TH1F("Eelectron_minus_ERecoJet", "E_{e}^{True} - E^{Reco}_{Jet} (|#eta^{Jet}|<0.7)",100,-20,30);
  //Detector Coordinate Histos
  TH1F * dPhiTj = new TH1F("dPhi_e_TrueJet", "|#Delta#varphi| (#varphi_{e} - #varphi^{True}_{Jet}) ", 32,0,M_PI);
  TH1F * dPhiRj = new TH1F("dPhi_e_RecoJet", "|#Delta#varphi| #varphi_{e} - #varphi(Jet^{Reco}_{Jet}) ", 32,0,M_PI);
  TH1F * dEtaTj = new TH1F("dEta_e_TrueJet", "|#Delta#eta| (#eta_{e} - #eta^{True}_{Jet})", 80,-10,10);
  TH1F * dEtaRj = new TH1F("dEta_e_RecoJet", "|#Delta#eta| (#eta_{e} - #eta^{Reco}_{Jet})", 80,-10,10);

  
  TH1F* Q2 = new TH1F("Q2","Q^{2}",500,0,500);
  TH1I * PID = new TH1I("PID_Histo", "PID",1000,-500,500);
  
  Long64_t nentries = _tree_event->GetEntries();
  for(Long64_t ie = 0; ie < nentries ; ie++){
    _tree_event->GetEntry(ie); //each entry contains a 5GeV Electron
    //fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ie, nentries);

    //Check for Truth Jets an NComponent > 2
    float Emin = 3.0;
    float hardest_jet_E = 0;
    int hardest = -1; //hardest jet index
    for (int j = 0; j < njets; j++) //ALL truth jet loop
      {
	if(isnan(e[j])) continue;
	justE->Fill(e[j]);
	if(matched_truthNComponent[j] < 2) continue; //cut on truth due to tracking lack of efficiency study
	if (e[j] < Emin) continue;
	
	if (e[j] > hardest_jet_E)
	    hardest_jet_E = e[j]; hardest = j;
	    
	std::cout<<"Jet # "<<j<<", number of Constituents = "<<matched_truthNComponent[j]<<std::endl;
	if (isnan(matched_truthE[j])) continue;
	  {
	    for (int c_index = 0; c_index < matched_truthNComponent[j]; c_index++) {
	      PID->Fill(all_Constituent_truthPID[j][c_index]);
	      // std::cout<<"all Constituent # "<<c_index<<" PID = "<<all_Constituent_truthPID[j][c_index];
	      // std::cout<<", All_Eta = "<<all_Constituent_truthEta[j][c_index];
	      // std::cout<<", Phi = "<<all_Constituent_truthPhi[j][c_index];
	      // std::cout<<", Pt = "<<all_Constituent_truthPt[j][c_index];
	      // std::cout<<", E = "<<all_Constituent_truthE[j][c_index]<<std::endl;
	    }
	    //std::cout<<std::endl;
	  }
      }
    //Detector Coordinate Histos
    Float_t True_DeltaPhi = TMath::Abs(TVector2::Phi_mpi_pi(electron_truthPhi - matched_truthPhi[hardest]));
    Float_t Reco_DeltaPhi = TMath::Abs(TVector2::Phi_mpi_pi(electron_truthPhi - phi[hardest]));
    dPhiTj->Fill(True_DeltaPhi);
    dPhiRj->Fill(Reco_DeltaPhi);
    dEtaTj->Fill(electron_truthEta-matched_truthEta[hardest]);
    dEtaRj->Fill(electron_truthEta-eta[hardest]);

    TLorentzVector e_vector;
    e_vector.SetPtEtaPhiE(electron_truthPt,electron_truthEta,electron_truthPhi,electron_truthE);
    Float_t Q_square = calc_Q_square(20,e_vector); //Electron Beam of 20 GeV/c
    Q2->Fill(Q_square);    
    //Inclusive Spectra
    Rjve->Fill(electron_truthE,e[hardest]);
    Tjve->Fill(electron_truthE,matched_truthE[hardest]);
    RjoTj->Fill(e[hardest]/matched_truthE[hardest]);
    
    //Kinematic Cuts
    // if (True_DeltaPhi < M_PI/2) continue;
    // if (matched_truthE[hardest] < 3.0) continue;
    // if (TMath::Abs(eta[hardest]) < 0.7) continue;
    //Electron/Jet Comparisons
    eoTj->Fill(matched_truthE[hardest]/electron_truthE);
    eoRj->Fill(e[hardest]/electron_truthE);
    emTj->Fill(electron_truthE-matched_truthE[hardest]);
    emRj->Fill(electron_truthE-e[hardest]);	
  } //entry loop

//Write to new root file
  TFile* fout = new TFile("Histograms_Jet_Callibration.root","RECREATE");
  justE->Write();
  PID->Write();
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

  RjoTj->GetXaxis()->SetTitle("E_{Jet}^{Reco}/E_{Jet}^{True}");
  RjoTj->Write();

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
  // H_dR->Write();
  // H_NExtra_Matches->Write();

  fout->Close();
  
}

