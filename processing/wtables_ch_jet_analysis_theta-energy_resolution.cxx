#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <filesystem>

// Root includes
#include <TROOT.h>
#include "TRint.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TVectorT.h"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

namespace fs = std::filesystem;
using namespace std;

// Forward-declaring functions
void prettyTH1F( TH1F * h1 , int color , int marker , float min , float max );
void better_yaxis(TH1F ** h1_array,int array_size);
TF1 * double_gaus(TH1F *h1,float min1, float max1, float min2, float max2,TString type ,int et, int p);
float constituent_pT_threshold(float eta);
TGraphErrors ** Th1_to_TGraph(TH1F ** h1_array, int array_size);
// ============================================================================================================================================
int main(int argc, char ** argv) {

#ifdef WITHRINT
  TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
#else
  TApplication *myapp = new TApplication("myapp",0,0);#
#endif

    if(argc < 5){
      cout << "Run this code as:\n\033[1;32m./analysis_momentum_resolution A B C D filename.root\033[0m\n";
      cout << "where:\nA = 1 -> Widths from table will be used\n  = 2 -> Widths from table \033[1;31mwon't\033[0m be used\n";
      cout << "B = 1 -> Table will be updated\n  = 2 -> Table \033[1;31mwon't\033[0m be updated\n";
      cout << "C = 1 -> Run code and quit\n  = 2 -> Run code and show plots\n";
      cout << "D = 1 -> Apply N_Missing < 1\n  = 2 -> Don't apply N_Missing Selection\n";
      exit(0);
    }

  bool use_widths = true;
  bool update_tab = false;
  bool keep_plots = false;
  bool N_Missing_Cut = true;
  int N_Missing_Max = 1;
  bool plot_mrad = false;
  string rad_mrad_string = "[rad]";
  if (plot_mrad)
    rad_mrad_string = "[mrad]";

  if (!fs::is_directory("../data"  ) || !fs::exists("../data"  )) fs::create_directory("../data"  ); // Create directory if it does not exist
  if (!fs::is_directory("tables"   ) || !fs::exists("tables"   )) fs::create_directory("tables"   );
  if (!fs::is_directory("../output") || !fs::exists("../output")) fs::create_directory("../output");
  if (!fs::is_directory("fits"     ) || !fs::exists("fits"     )) fs::create_directory("fits"     );
  if (!fs::is_directory("results"  ) || !fs::exists("results"  )) fs::create_directory("results"  );

  cout << "\033[1;31m********************************************************************\nUSEFUL INFO:\033[0m\nWill be loading data from file: '" << argv[5] << "' assumed to be in directory 'data'" << endl;

  if     (atoi(argv[1])==1){use_widths = true ;	cout << "Will be using widths from table\n" ;}
  else if(atoi(argv[1])==2){use_widths = false;	cout << "Won't be using widths from table\n";}
  else{cout << "Something wrong with your election of input parameter 'A'. Bailing out!\n"; exit(0);}

  if     (atoi(argv[2])==1){update_tab = true ;   cout << "Table will be updated\n" ;}
  else if(atoi(argv[2])==2){update_tab = false;   cout << "Table won't be updated\n";}
  else{cout << "Something wrong with your election of input parameter 'B'. Bailing out!\n"; exit(0);}

  if     (atoi(argv[3])==1){keep_plots = false;	cout << "Will run and quit. Examine the output files for resulting plots\n";gROOT->SetBatch(kTRUE);}
  else if(atoi(argv[3])==2){keep_plots = true ;	cout << "Will run and show the plots\n" ; }
  else{cout << "Something wrong with your election of input parameter 'C'. Bailing out!\n"; exit(0);}

  if     (atoi(argv[4])==1){N_Missing_Cut = true ;   cout << "Will apply N Missing Const. Cut\n" ;}
  else if(atoi(argv[4])==2){N_Missing_Cut = false;   cout << "Will NOT apply N Missing Const. Cut\n";}
  else{cout << "Something wrong with your election of input parameter 'D'. Bailing out!\n"; exit(0);}


  // -------------------------
  // Binning
  float eta_bin[] = {-3.0,-1.5,-.5,0.5,1.5,3.0};
  //float eta_bin[] = {-3.0,-1.0,0.0,1.0,3.0};
  float mom_bin[] = {4.,6.,8.,10.,12.,15, 20.};

  const int size_eta_bin = sizeof(eta_bin)/sizeof(*eta_bin);
  const int size_mom_bin = sizeof(mom_bin)/sizeof(*mom_bin);

  TVectorT<double> TVT_eta_bin(size_eta_bin);	for(int i = 0 ; i < size_eta_bin ; i++) TVT_eta_bin[i] = eta_bin[i];
  TVectorT<double> TVT_mom_bin(size_mom_bin);	for(int i = 0 ; i < size_mom_bin ; i++) TVT_mom_bin[i] = mom_bin[i];
  // -------------------------
  // useful strings
  string raw_fname = argv[5];
  TString infile = "../data/" + raw_fname;
  raw_fname.resize(raw_fname.size()-5);
  TString outfile = "../output/No_Missing_const_output_mom_res_" + raw_fname + Form("sigma_eta_%i_p_%i_",size_eta_bin-1,size_mom_bin-1) + ".root";
  if (!N_Missing_Cut)
    outfile = "../output/Missing_const_output_mom_res_" + raw_fname + Form("sigma_eta_%i_p_%i_",size_eta_bin-1,size_mom_bin-1) + ".root";
  TString tab_name = "tables/tab_mom_res_" + raw_fname + Form("sigma_eta_%i_p_%i_",size_eta_bin-1,size_mom_bin-1) + ".txt";
  TString out_pdf = "output_fits_mom_res_" + raw_fname + Form("sigma_eta_%i_p_%i_",size_eta_bin-1,size_mom_bin-1) + ".pdf";
  TString out_pdf2 = "results/results_mom_res_" + raw_fname + Form("sigma_eta_%i_p_%i_",size_eta_bin-1,size_mom_bin-1) + ".pdf";
  // -------------------------------------------------------------
  // Some settings
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  gStyle -> SetOptStat(0);	
  // -------------------------------------------------------------
  // Loading all the needed info from the root file
  TFile * F = new TFile(infile);
  TTree *T = dynamic_cast<TTree *>(F->Get("T"));
  if (T == NULL) { std::cout << " Tree Fail " << std::endl; exit(EXIT_FAILURE); }
  Int_t njets;
  int nEntries = T -> GetEntries();
  const int MaxNumJets = 20;
  const int kMaxConstituents = 100;
  array<Float_t, MaxNumJets> E,Eta,Phi,Pt,gE,gEta,gPhi,gPt;
  array<Int_t, MaxNumJets> NComponent, gNComponent;
  Float_t electron_gE,electron_gEta,electron_gPhi,electron_gPt;
  array<array<Float_t, kMaxConstituents >, MaxNumJets > gComponent_Eta,gComponent_PID,
    gComponent_Pt,gComponent_Phi,gComponent_E, gComponent_Charge;

  T -> SetBranchAddress("njets",&njets);
  T -> SetBranchAddress("e",&E);
  T -> SetBranchAddress("eta",&Eta);
  T -> SetBranchAddress("phi",&Phi);
  T -> SetBranchAddress("pt",&Pt);
  T -> SetBranchAddress("nComponent",&NComponent);

  T -> SetBranchAddress("matched_truthE",&gE);
  T -> SetBranchAddress("matched_truthEta",&gEta);
  T -> SetBranchAddress("matched_truthPhi",&gPhi);
  T -> SetBranchAddress("matched_truthPt",&gPt);
  T -> SetBranchAddress("matched_truthNComponent",&gNComponent);

  T -> SetBranchAddress("matched_Constituent_truthEta", gComponent_Eta.data());
  T -> SetBranchAddress("matched_Constituent_truthPID",gComponent_PID.data());
  T -> SetBranchAddress("matched_Constituent_truthPt",gComponent_Pt.data());
  T -> SetBranchAddress("matched_Constituent_truthE",gComponent_E.data());
  T -> SetBranchAddress("matched_Constituent_truthPhi",gComponent_Phi.data());
  T -> SetBranchAddress("matched_Constituent_truthCharge",gComponent_Charge.data());

  // -------------------------------------------------------------
  fstream tab;
  
  float approx_sig_dpp[size_eta_bin-1][size_mom_bin-1] = {{0}};
  float approx_sig_dth[size_eta_bin-1][size_mom_bin-1] = {{0}};
  float approx_sig_dph[size_eta_bin-1][size_mom_bin-1] = {{0}};
  float approx_sig_des[size_eta_bin-1][size_mom_bin-1] = {{0}};

  float approx_mean_dpp[size_eta_bin-1][size_mom_bin-1] = {{0}};
  float approx_mean_dth[size_eta_bin-1][size_mom_bin-1] = {{0}};
  float approx_mean_dph[size_eta_bin-1][size_mom_bin-1] = {{0}};
  float approx_mean_des[size_eta_bin-1][size_mom_bin-1] = {{0}};

  TString temp_str;
  if(use_widths){
    tab.open(tab_name);
    if(!tab){cout << "Could not find file '" << tab_name << "'" << endl; use_widths = false; update_tab = true;}
    else{
      cout << "Loading parameters from file '" << tab_name << "'" << endl;
      for(int et = 0 ; et < size_eta_bin-1 ; et++){ for(int p = 0 ; p < size_mom_bin-1 ; p++){tab >> approx_sig_dpp[et][p];}}
      for(int et = 0 ; et < size_eta_bin-1 ; et++){ for(int p = 0 ; p < size_mom_bin-1 ; p++){tab >> approx_sig_dth[et][p];}}
      for(int et = 0 ; et < size_eta_bin-1 ; et++){ for(int p = 0 ; p < size_mom_bin-1 ; p++){tab >> approx_sig_dph[et][p];}}
      for(int et = 0 ; et < size_eta_bin-1 ; et++){ for(int p = 0 ; p < size_mom_bin-1 ; p++){tab >> approx_sig_des[et][p];}}
      for(int et = 0 ; et < size_eta_bin-1 ; et++){ for(int p = 0 ; p < size_mom_bin-1 ; p++){tab >> approx_mean_dpp[et][p];}}
      for(int et = 0 ; et < size_eta_bin-1 ; et++){ for(int p = 0 ; p < size_mom_bin-1 ; p++){tab >> approx_mean_dth[et][p];}}
      for(int et = 0 ; et < size_eta_bin-1 ; et++){ for(int p = 0 ; p < size_mom_bin-1 ; p++){tab >> approx_mean_dph[et][p];}}
      for(int et = 0 ; et < size_eta_bin-1 ; et++){ for(int p = 0 ; p < size_mom_bin-1 ; p++){tab >> approx_mean_des[et][p];}}
    }
    tab.close();
  }

  float approx_sig_dpp_3_0[size_eta_bin-1][size_mom_bin-1] = {{0}}; float approx_sig_dpp_1_1[size_eta_bin-1][size_mom_bin-1] = {{0}};
  float approx_sig_dth_3_0[size_eta_bin-1][size_mom_bin-1] = {{0}}; float approx_sig_dth_1_1[size_eta_bin-1][size_mom_bin-1] = {{0}};
  float approx_sig_dph_3_0[size_eta_bin-1][size_mom_bin-1] = {{0}}; float approx_sig_dph_1_1[size_eta_bin-1][size_mom_bin-1] = {{0}};
  float approx_sig_des_3_0[size_eta_bin-1][size_mom_bin-1] = {{0}}; float approx_sig_des_1_1[size_eta_bin-1][size_mom_bin-1] = {{0}};

  for(int et = 0 ; et < size_eta_bin-1 ; et++){
    for(int p = 0 ; p < size_mom_bin-1 ; p++){
      approx_sig_dpp_3_0[et][p] = 3.0*approx_sig_dpp[et][p]; approx_sig_dpp_1_1[et][p] = 1.1*approx_sig_dpp[et][p];
      approx_sig_dth_3_0[et][p] = 3.0*approx_sig_dth[et][p]; approx_sig_dth_1_1[et][p] = 1.1*approx_sig_dth[et][p];
      approx_sig_dph_3_0[et][p] = 3.0*approx_sig_dph[et][p]; approx_sig_dph_1_1[et][p] = 1.1*approx_sig_dph[et][p];
      approx_sig_des_3_0[et][p] = 3.0*approx_sig_des[et][p]; approx_sig_des_1_1[et][p] = 1.1*approx_sig_des[et][p];
    }
  }
  // -------------------------------------------------------------
  // Defining histograms
  TH1F *** h1_dpp_p_et_bins = new TH1F**[size_eta_bin-1];	// delta p / p vs. p in eta bins
  TH1F *** h1_dth_p_et_bins = new TH1F**[size_eta_bin-1];	// delta theta vs. p in eta bins
  TH1F *** h1_dph_p_et_bins = new TH1F**[size_eta_bin-1];	// delta phi   vs. p in eta bins
  TH1F *** h1_des_p_et_bins = new TH1F**[size_eta_bin-1];	// delta phi   vs. p in eta bins

  int n_dpp_bins = 30;
  int n_dth_bins = 30;
  int n_dph_bins = 30;
  int n_des_bins = 30;

  for(int et = 0 ; et < size_eta_bin-1 ; et++){
    h1_dpp_p_et_bins[et] = new TH1F*[size_mom_bin-1];
    h1_dth_p_et_bins[et] = new TH1F*[size_mom_bin-1];
    h1_dph_p_et_bins[et] = new TH1F*[size_mom_bin-1];
    h1_des_p_et_bins[et] = new TH1F*[size_mom_bin-1];

    for(int p = 0 ; p < size_mom_bin-1 ; p++){

      if(use_widths){
        h1_dpp_p_et_bins[et][p] = new TH1F(Form("h1_dpp_p_et_bins_%i_%i",et,p),";dP/P;Counts",n_dpp_bins,
            approx_mean_dpp[et][p]-approx_sig_dpp_3_0[et][p],approx_mean_dpp[et][p]+approx_sig_dpp_3_0[et][p]);

        h1_dth_p_et_bins[et][p] = new TH1F(Form("h1_dth_p_et_bins_%i_%i",et,p),";d#theta [rad];Counts",n_dth_bins,
            approx_mean_dth[et][p]-approx_sig_dth_3_0[et][p],approx_mean_dth[et][p]+approx_sig_dth_3_0[et][p]);

        h1_dph_p_et_bins[et][p] = new TH1F(Form("h1_dph_p_et_bins_%i_%i",et,p),";d#phi [rad];Counts",n_dph_bins,
            approx_mean_dph[et][p]-approx_sig_dph_3_0[et][p],approx_mean_dph[et][p]+approx_sig_dph_3_0[et][p]);

        h1_des_p_et_bins[et][p] = new TH1F(Form("h1_des_p_et_bins_%i_%i",et,p),";d(E_{Reco}/E_{Truth});Counts",n_des_bins,
            approx_mean_des[et][p]-approx_sig_des_3_0[et][p],approx_mean_des[et][p]+approx_sig_des_3_0[et][p]);
      }

      else if (N_Missing_Cut){
        h1_dpp_p_et_bins[et][p] = new TH1F(Form("h1_dpp_p_et_bins_%i_%i",et,p),
            ";dP/P;Counts",n_dpp_bins,-0.06,0.06);
        h1_dth_p_et_bins[et][p] = new TH1F(Form("h1_dth_p_et_bins_%i_%i",et,p),
            ";d#theta [rad];Counts",n_dth_bins,-0.005,0.005);
        h1_dph_p_et_bins[et][p] = new TH1F(Form("h1_dph_p_et_bins_%i_%i",et,p),
            ";d#phi [rad];Counts"  ,n_dph_bins,-0.01,0.01);
        h1_des_p_et_bins[et][p] = new TH1F(Form("h1_des_p_et_bins_%i_%i",et,p),
            ";Energy Scale;Counts",n_des_bins,0.75,1.05);
      }

      else {
        h1_dpp_p_et_bins[et][p] = new TH1F(Form("h1_dpp_p_et_bins_%i_%i",et,p),
            ";dP/P;Counts",60,-0.4,0.8);
        h1_dth_p_et_bins[et][p] = new TH1F(Form("h1_dth_p_et_bins_%i_%i",et,p),
            ";d#theta [rad];Counts",n_dth_bins,-0.2,0.2);
        h1_dph_p_et_bins[et][p] = new TH1F(Form("h1_dph_p_et_bins_%i_%i",et,p),
            ";d#phi [rad];Counts",n_dph_bins,-0.2,0.2);
        h1_des_p_et_bins[et][p] = new TH1F(Form("h1_des_p_et_bins_%i_%i",et,p),
            ";Energy Scale;Counts",n_des_bins,0.0,2.);
      }

      h1_dpp_p_et_bins[et][p] -> SetTitle(Form("%.1f < |#eta| < %.1f, %.1f < P < %.1f GeV/c",eta_bin[et],eta_bin[et+1],mom_bin[p],mom_bin[p+1]));
      h1_dth_p_et_bins[et][p] -> SetTitle(Form("%.1f < |#eta| < %.1f, %.1f < P < %.1f GeV/c",eta_bin[et],eta_bin[et+1],mom_bin[p],mom_bin[p+1]));
      h1_dph_p_et_bins[et][p] -> SetTitle(Form("%.1f < |#eta| < %.1f, %.1f < P < %.1f GeV/c",eta_bin[et],eta_bin[et+1],mom_bin[p],mom_bin[p+1]));
      h1_des_p_et_bins[et][p] -> SetTitle(Form("%.1f < |#eta| < %.1f, %.1f < E < %.1f GeV",eta_bin[et],eta_bin[et+1],mom_bin[p],mom_bin[p+1]));
    } 
  }	
  // -------------------------------------------------------------	
  TH1F ** h1_dpp_v_p_et_bins = new TH1F*[size_eta_bin-1];
  TH1F ** h1_dth_v_p_et_bins = new TH1F*[size_eta_bin-1];
  TH1F ** h1_dph_v_p_et_bins = new TH1F*[size_eta_bin-1];
  TH1F ** h1_des_v_p_et_bins = new TH1F*[size_eta_bin-1];
  TH2F * mom_response = new TH2F("mom_response","P_{Charge}^{Reco} vs. P_{Charge}^{Truth}",50,0,50,50,0,50);

  for(int et = 0 ; et < size_eta_bin-1 ; et++){
    if (N_Missing_Cut) h1_dpp_v_p_et_bins[et] = new TH1F(Form("h1_dpp_v_p_et_bins_%i",et),";P [GeV/c];dP/P [%]",size_mom_bin-1,mom_bin);
    else h1_dpp_v_p_et_bins[et] = new TH1F(Form("h1_dpp_v_p_et_bins_%i",et),";P [GeV/c];StdDev dP/P [%]",size_mom_bin-1,mom_bin);
    prettyTH1F( h1_dpp_v_p_et_bins[et] , et , 20 , 1. , 100. );

    h1_dth_v_p_et_bins[et] = new TH1F(Form("h1_dth_v_p_et_bins_%i",et),Form(";P [GeV/c];d#theta %s",rad_mrad_string.c_str()),size_mom_bin-1,mom_bin);	prettyTH1F( h1_dth_v_p_et_bins[et] , et , 20 , 0.001 , 1.  );
    h1_dph_v_p_et_bins[et] = new TH1F(Form("h1_dph_v_p_et_bins_%i",et),Form(";P [GeV/c];d#phi %s",rad_mrad_string.c_str())  ,size_mom_bin-1,mom_bin);	prettyTH1F( h1_dph_v_p_et_bins[et] , et , 20 , 0.001 , 1. );
    h1_des_v_p_et_bins[et] = new TH1F(Form("h1_des_v_p_et_bins_%i",et),";E [GeV];d(E_{Reco}/E_{True})",size_mom_bin-1,mom_bin);	//prettyTH1F( h1_des_v_p_et_bins[et] , et , 20 , 1., 100. );
  }

  TH1F ** h1_dpp_v_et_p_bins = new TH1F*[size_mom_bin-1];
  TH1F ** h1_dth_v_et_p_bins = new TH1F*[size_mom_bin-1];
  TH1F ** h1_dph_v_et_p_bins = new TH1F*[size_mom_bin-1];
  TH1F ** h1_des_v_et_p_bins = new TH1F*[size_mom_bin-1];

  for(int p = 0 ; p < size_mom_bin-1 ; p++){
    if (N_Missing_Cut) h1_dpp_v_et_p_bins[p] = new TH1F(Form("h1_dpp_v_et_p_bins_%i",p),";#eta;dP/P [%]"      ,size_eta_bin-1,eta_bin);
    else h1_dpp_v_et_p_bins[p] = new TH1F(Form("h1_dpp_v_et_p_bins_%i",p),";P [GeV/c];StdDev dP/P [%]",size_eta_bin-1,mom_bin);
    prettyTH1F( h1_dpp_v_et_p_bins[p] , p , 20 , 1. , 100. );
    
    h1_dth_v_et_p_bins[p] = new TH1F(Form("h1_dth_v_et_p_bins_%i",p),Form(";#eta;d#theta %s",rad_mrad_string.c_str()),size_eta_bin-1,eta_bin);
    prettyTH1F( h1_dth_v_et_p_bins[p] , p , 20 , 0.001 , 1.  );

    h1_dph_v_et_p_bins[p] = new TH1F(Form("h1_dph_v_et_p_bins_%i",p),Form(";#eta;d#phi [rad]",rad_mrad_string.c_str()),size_eta_bin-1,eta_bin);
    prettyTH1F( h1_dph_v_et_p_bins[p] , p , 20 , 0.001 , 1000. );

    h1_des_v_et_p_bins[p] = new TH1F(Form("h1_des_v_et_p_bins_%i",p),";#eta;d(E_{Reco}/E_{True})",size_eta_bin-1,eta_bin);
    prettyTH1F( h1_des_v_et_p_bins[p] , p , 20 , 1. , 100. );
  }

  cout << "\033[1;31m********************************************************************\033[0m\n";
  // -------------------------------------------------------------
  // Loop over entries of the tree	
  for (Long64_t ev = 0; ev < nEntries; ev++){
    T->GetEntry(ev);
    if (ev%10000==0) fprintf(stderr,"\r%d: Entry %lli out of %d",__LINE__,ev,nEntries);
    //if (ev==500000) break;
    for (int n = 0; n < njets; ++n) {

      if (NComponent[n] < 4) continue;
      if (isnan(gE[n])) continue;
      ROOT::Math::PtEtaPhiEVector Lorentz(Pt[n],Eta[n],Phi[n],E[n]);
      ROOT::Math::PtEtaPhiEVector gLorentz(gPt[n],gEta[n],gPhi[n],gE[n]);

      bool eta_const_cut = true; //avoid crack where barrel meets disks
      bool pt_const_cut = true;//cut away helixes
      int n_neutral = 0;
      int n_ch = 0;
      int n_neutral_max = 1;
      bool maxed_neutrals = true;
      for (int i = 0; i < gNComponent[n]; i++){

        eta_const_cut = (  ( (abs(gComponent_Eta[n][i]) > 1.06) && (abs(gComponent_Eta[n][i]) < 1.13) )
            || (abs(gComponent_Eta[n][i]) > 3.5)  );

        pt_const_cut = (gComponent_Pt[n][i] < constituent_pT_threshold(gComponent_Eta[n][i]));

        if (eta_const_cut || pt_const_cut) break; //skip jets that fail (general cut)

        if (gComponent_Charge[n][i] == 0)
          n_neutral++;
        else
          n_ch++;

        ROOT::Math::PtEtaPhiEVector gConstLorentz(gComponent_Pt[n][i],
            gComponent_Eta[n][i],
            gComponent_Phi[n][i],
            gComponent_E[n][i]);
        if (gComponent_Charge[n][i] == 0)
          gLorentz -= gConstLorentz;
      }

      int N_Missing = gNComponent[n] - NComponent[n] - n_neutral;
      if ( (N_Missing_Cut) && (N_Missing >= N_Missing_Max) ) continue;
      if ( (!N_Missing_Cut) && (N_Missing < N_Missing_Max) ) continue; //anticut of above
      float dth = Lorentz.Theta() - gLorentz.Theta();
      float geta = gLorentz.Eta();
      float P_reco = Lorentz.P();
      float P_truth = gLorentz.P();
      float dP_P = (P_truth-P_reco)/P_truth;

      float E_gE = (Lorentz.E()/gLorentz.E());
      //float E_gE = (E[n]/gE[n]);
      mom_response->Fill(gLorentz.P(),Lorentz.P());

      float dph = Lorentz.Phi() - gLorentz.Phi();
      float dR = ROOT::Math::VectorUtil::DeltaR(Lorentz,gLorentz);

      //Filling 
      for(int et = 0 ; et < size_eta_bin-1 ; et++){
        if( geta >  eta_bin[et] &&  geta <= eta_bin[et+1] ){
          for(int p = 0 ; p < size_mom_bin-1 ; p++){
            if( P_truth > mom_bin[p] && P_truth <= mom_bin[p+1] ){
              h1_dpp_p_et_bins[et][p] -> Fill( dP_P );
              h1_dth_p_et_bins[et][p] -> Fill( dth  );
              h1_dph_p_et_bins[et][p] -> Fill( dph  );
            }	
            if ((gLorentz.E() > mom_bin[p]) && (gLorentz.E() <= mom_bin[p+1]))
              h1_des_p_et_bins[et][p] -> Fill( E_gE );
            //if (gE[n] > mom_bin[p] && gE[n] <= mom_bin[p+1])//P != E, but binning is close enough:w
            //h1_des_p_et_bins[et][p] -> Fill( E_gE );
          }
        }
      }
    }
  }
  cout << "\033[1;31m********************************************************************\033[0m\n";
  // -------------------------------------------------------------
  // Declaring other useful variables and functions
  float width_dpp[size_eta_bin-1][size_mom_bin-1] = {{0}};
  float error_dpp[size_eta_bin-1][size_mom_bin-1] = {{0}};

  float width_dth[size_eta_bin-1][size_mom_bin-1] = {{0}};
  float error_dth[size_eta_bin-1][size_mom_bin-1] = {{0}};

  float width_dph[size_eta_bin-1][size_mom_bin-1] = {{0}};
  float error_dph[size_eta_bin-1][size_mom_bin-1] = {{0}};

  float width_des[size_eta_bin-1][size_mom_bin-1] = {{0}};
  float error_des[size_eta_bin-1][size_mom_bin-1] = {{0}};

  float mean_dpp[size_eta_bin-1][size_mom_bin-1] = {{0}};
  float mean_dth[size_eta_bin-1][size_mom_bin-1] = {{0}};
  float mean_dph[size_eta_bin-1][size_mom_bin-1] = {{0}};
  float mean_des[size_eta_bin-1][size_mom_bin-1] = {{0}};

  TF1 *** f_gaus_dpp = new TF1**[size_eta_bin-1];
  TF1 *** f_gaus_dth = new TF1**[size_eta_bin-1];
  TF1 *** f_gaus_dph = new TF1**[size_eta_bin-1];
  TF1 *** f_gaus_des = new TF1**[size_eta_bin-1];

  for(int et = 0 ; et < size_eta_bin-1 ; et++){
    f_gaus_dpp[et] = new TF1*[size_mom_bin-1];
    f_gaus_dth[et] = new TF1*[size_mom_bin-1];
    f_gaus_dph[et] = new TF1*[size_mom_bin-1];
    f_gaus_des[et] = new TF1*[size_mom_bin-1];

    for(int p = 0 ; p < size_mom_bin-1 ; p++){
      if(use_widths){

        f_gaus_dpp[et][p] = new TF1(Form("f_gaus_dpp_%i_%i",et,p),"gaus",approx_mean_dpp[et][p]-approx_sig_dpp_1_1[et][p],
            approx_mean_dpp[et][p]+approx_sig_dpp_1_1[et][p]);
        f_gaus_dth[et][p] = new TF1(Form("f_gaus_dth_%i_%i",et,p),"gaus",approx_mean_dth[et][p]-approx_sig_dth_1_1[et][p],
            approx_mean_dth[et][p]+approx_sig_dth_1_1[et][p]);
        f_gaus_dph[et][p] = new TF1(Form("f_gaus_dph_%i_%i",et,p),"gaus",approx_mean_dph[et][p]-approx_sig_dph_1_1[et][p],
            approx_mean_dph[et][p]+approx_sig_dph_1_1[et][p]);
        f_gaus_des[et][p] = new TF1(Form("f_gaus_des_%i_%i",et,p),"gaus",approx_mean_des[et][p]-approx_sig_des_1_1[et][p],
            approx_mean_des[et][p]+approx_sig_des_1_1[et][p]);
      }

      else if (N_Missing_Cut){
          
        f_gaus_dpp[et][p] = new TF1(Form("f_gaus_dpp_%i_%i",et,p),"gaus",-0.01,0.01);
        f_gaus_dth[et][p] = new TF1(Form("f_gaus_dth_%i_%i",et,p),"gaus",-0.0008,0.0008);
        f_gaus_dph[et][p] = new TF1(Form("f_gaus_dph_%i_%i",et,p),"gaus",-0.002,0.002);
        f_gaus_des[et][p] = new TF1(Form("f_gaus_des_%i_%i",et,p),"gaus",0.5,1.0);
      }

      else{//n_miss
        f_gaus_dpp[et][p] = new TF1(Form("f_gaus_dpp_%i_%i",et,p),"gaus",-0.02,0.02);
        f_gaus_dth[et][p] = new TF1(Form("f_gaus_dth_%i_%i",et,p),"gaus",-0.03,0.03);
        f_gaus_dph[et][p] = new TF1(Form("f_gaus_dph_%i_%i",et,p),"gaus",-0.03,0.03);
        f_gaus_des[et][p] = new TF1(Form("f_gaus_des_%i_%i",et,p),"gaus",0.5,1.0);
      }
    }//p
  }//et

  // -------------------------------------------------------------
  // Doing fits
  TCanvas ** c_fits_p  = new TCanvas*[size_eta_bin-1];
  TCanvas ** c_fits_th = new TCanvas*[size_eta_bin-1];
  TCanvas ** c_fits_ph = new TCanvas*[size_eta_bin-1];
  TCanvas ** c_fits_es = new TCanvas*[size_eta_bin-1];

  for(int et = 0 ; et < size_eta_bin-1 ; et++){
    c_fits_p [et] = new TCanvas(Form("c_fits_p_%i" ,et),Form("dP/P  , %.1f<eta<%.1f",eta_bin[et],eta_bin[et+1]),1000,800);
    c_fits_p [et] -> Divide(5,2);
    c_fits_th[et] = new TCanvas(Form("c_fits_th_%i",et),Form("dtheta, %.1f<eta<%.1f",eta_bin[et],eta_bin[et+1]),1000,800);
    c_fits_th[et] -> Divide(5,2);
    c_fits_ph[et] = new TCanvas(Form("c_fits_ph_%i",et),Form("dphi  , %.1f<eta<%.1f",eta_bin[et],eta_bin[et+1]),1000,800);
    c_fits_ph[et] -> Divide(5,2);
    c_fits_es[et] = new TCanvas(Form("c_fits_es_%i",et),Form("d(E_{Reco}/E_{Truth}), %.1f<eta<%.1f",eta_bin[et],eta_bin[et+1]),1000,800);
    c_fits_es[et] -> Divide(5,2);

    for(int p = 0 ; p < size_mom_bin-1 ; p++){

      // Momentum resolutions
      c_fits_p [et] -> cd(p+1);
      h1_dpp_p_et_bins[et][p] -> Draw();	//h1_dpp_p_et_bins[et][p] -> Fit(Form("f_gaus_dpp_%i_%i",et,p),"RQ");

      if(N_Missing_Cut){
        h1_dpp_p_et_bins[et][p] -> Fit(f_gaus_dpp[et][p],"RQ");			
        width_dpp[et][p] = f_gaus_dpp[et][p] -> GetParameter(2);
        error_dpp[et][p] = (f_gaus_dpp[et][p] -> GetParError(2));
        mean_dpp[et][p] = f_gaus_dpp[et][p] -> GetParameter(1);
      }

      else{
        width_dpp[et][p] = h1_dpp_p_et_bins[et][p]->GetStdDev();
        error_dpp[et][p] = width_dpp[et][p]/sqrt(2*(h1_dpp_p_et_bins[et][p]->GetEntries())-2);
        mean_dpp[et][p] = h1_dpp_p_et_bins[et][p]->GetMean();

        //Full Width at Half Max for dP/P
        //double max = h1_dpp_p_et_bins[et][p]->GetMaximum();
        //double fwhm_left = h1_dpp_p_et_bins[et][p]->FindFirstBinAbove(max/2.);
        //double fwhm_right = h1_dpp_p_et_bins[et][p]->FindLastBinAbove(max/2.);
        //double fwhm = h1_dpp_p_et_bins[et][p]->GetBinCenter(fwhm_right) - h1_dpp_p_et_bins[et][p]->GetBinCenter(fwhm_left); 
        //width_dpp[et][p] = fwhm;
        //cout<<"max = "<<max<<" left = "<<fwhm_left<<" right = "<<fwhm_right<<" width = "<<width_dpp[et][p]<<endl;
      }

      // Theta resolution
      c_fits_th[et] -> cd(p+1);
      float rad_mrad = 1.;;
      if (plot_mrad)
        rad_mrad = 1000.;
      h1_dth_p_et_bins[et][p] -> Draw();	h1_dth_p_et_bins[et][p] -> Fit(Form("f_gaus_dth_%i_%i",et,p),"RQ");
      width_dth[et][p] = f_gaus_dth[et][p] -> GetParameter(2)*rad_mrad;
      error_dth[et][p] = f_gaus_dth[et][p] -> GetParError(2)*rad_mrad;
      mean_dth[et][p] = f_gaus_dth[et][p] -> GetParameter(1)*rad_mrad;

      // Phi resolution
      c_fits_ph[et] -> cd(p+1);
      h1_dph_p_et_bins[et][p] -> Draw();	h1_dph_p_et_bins[et][p] -> Fit(Form("f_gaus_dph_%i_%i",et,p),"RQ");
      width_dph[et][p] = f_gaus_dph[et][p] -> GetParameter(2)*rad_mrad;
      error_dph[et][p] = (f_gaus_dph[et][p] -> GetParError(2)*rad_mrad);
      mean_dph[et][p] = f_gaus_dph[et][p] -> GetParameter(1);

      // Energy Scale resolution
      c_fits_es[et] -> cd(p+1);
      h1_des_p_et_bins[et][p] -> Draw();	h1_des_p_et_bins[et][p] -> Fit(Form("f_gaus_des_%i_%i",et,p),"RQ");
      width_des[et][p] = f_gaus_des[et][p] -> GetParameter(2);
      error_des[et][p] = (f_gaus_des[et][p] -> GetParError(2));
      mean_des[et][p] = f_gaus_des[et][p] -> GetParameter(1);

      cout<<__LINE__<<": dph width = "<<width_dph[et][p]<<endl;
      cout<<__LINE__<<": dpp width = "<<width_dpp[et][p]<<endl;
      cout<<__LINE__<<": dth width = "<<width_dth[et][p]<<endl;
      cout<<__LINE__<<": des width = "<<width_des[et][p]<<endl;
      // ----
      if(h1_dpp_p_et_bins[et][p]->GetMaximum()>50.){
        h1_dpp_v_p_et_bins[et] -> SetBinContent(p +1,width_dpp[et][p]*100. ); //*100
        h1_dpp_v_p_et_bins[et] -> SetBinError  (p +1,error_dpp[et][p]*100. ); 
        cout<<"Bin set to "<<h1_dpp_v_p_et_bins[et] -> GetBinContent(p +1)<<endl;;
      }
      if(h1_dth_p_et_bins[et][p]->GetMaximum()>50.){
        h1_dth_v_p_et_bins[et] -> SetBinContent(p +1,width_dth[et][p]);
        h1_dth_v_p_et_bins[et] -> SetBinError  (p +1,error_dth[et][p]);
      }
      if(h1_dph_p_et_bins[et][p]->GetMaximum()>50.){
        h1_dph_v_p_et_bins[et] -> SetBinContent(p +1,width_dph[et][p]);
        h1_dph_v_p_et_bins[et] -> SetBinError  (p +1,error_dph[et][p]);
      }
      if(h1_des_p_et_bins[et][p]->GetMaximum()>50.){
        h1_des_v_p_et_bins[et] -> SetBinContent(p +1,width_des[et][p]);
        h1_des_v_p_et_bins[et] -> SetBinError  (p +1,error_des[et][p]);
      }
      if(h1_dpp_p_et_bins[et][p]->GetMaximum()>50.){
        h1_dpp_v_et_p_bins[ p] -> SetBinContent(et+1,width_dpp[et][p]*100. );
        h1_dpp_v_et_p_bins[ p] -> SetBinError  (et+1,error_dpp[et][p]*100. );
      }
      if(h1_dth_p_et_bins[et][p]->GetMaximum()>50.){
        h1_dth_v_et_p_bins[ p] -> SetBinContent(et+1,width_dth[et][p]);
        h1_dth_v_et_p_bins[ p] -> SetBinError  (et+1,error_dth[et][p]);
      }
      if(h1_dph_p_et_bins[et][p]->GetMaximum()>50.){
        h1_dph_v_et_p_bins[ p] -> SetBinContent(et+1,width_dph[et][p]);
        h1_dph_v_et_p_bins[ p] -> SetBinError  (et+1,error_dph[et][p]);
      }
      if(h1_des_p_et_bins[et][p]->GetMaximum()>50.){
        h1_des_v_et_p_bins[ p] -> SetBinContent(et+1,width_des[et][p]);
        h1_des_v_et_p_bins[ p] -> SetBinError  (et+1,error_des[et][p]);
      }

    }

    c_fits_p [et] -> Modified();	c_fits_p [et] -> Update();
    c_fits_th[et] -> Modified();	c_fits_th[et] -> Update();
    c_fits_ph[et] -> Modified();	c_fits_ph[et] -> Update();
    c_fits_es[et] -> Modified();	c_fits_es[et] -> Update();
  }
  // -------------------------------------------------------------
  // Updating table with width values
  ofstream updated_tab;
  if(update_tab){
    updated_tab.open(tab_name);
    //widths
    for(int et = 0 ; et < size_eta_bin-1 ; et++){
      for(int p = 0 ; p < size_mom_bin-1 ; p++){
        updated_tab << width_dpp[et][p]; //what to do if the fit did not converge?
        if(p == size_mom_bin-2) updated_tab << "\n";
        else updated_tab << "\t";
      }
    }
    updated_tab << "\n";
    for(int et = 0 ; et < size_eta_bin-1 ; et++){
      for(int p = 0 ; p < size_mom_bin-1 ; p++){
        updated_tab << width_dth[et][p];
        if(p == size_mom_bin-2) updated_tab << "\n";
        else updated_tab << "\t";
      }
    }
    updated_tab << "\n";
    for(int et = 0 ; et < size_eta_bin-1 ; et++){
      for(int p = 0 ; p < size_mom_bin-1 ; p++){
        updated_tab << width_dph[et][p];
        if(p == size_mom_bin-2) updated_tab << "\n";
        else updated_tab << "\t";
      }
    }
    updated_tab << "\n";
    for(int et = 0 ; et < size_eta_bin-1 ; et++){
      for(int p = 0 ; p < size_mom_bin-1 ; p++){
        updated_tab << width_des[et][p];
        if(p == size_mom_bin-2) updated_tab << "\n";
        else updated_tab << "\t";
      }
    }
    //Means
    for(int et = 0 ; et < size_eta_bin-1 ; et++){
      for(int p = 0 ; p < size_mom_bin-1 ; p++){
        updated_tab << mean_dpp[et][p];
        if(p == size_mom_bin-2) updated_tab << "\n";
        else updated_tab << "\t";
      }
    }
    updated_tab << "\n";
    for(int et = 0 ; et < size_eta_bin-1 ; et++){
      for(int p = 0 ; p < size_mom_bin-1 ; p++){
        updated_tab << mean_dth[et][p];
        if(p == size_mom_bin-2) updated_tab << "\n";
        else updated_tab << "\t";
      }
    }
    updated_tab << "\n";
    for(int et = 0 ; et < size_eta_bin-1 ; et++){
      for(int p = 0 ; p < size_mom_bin-1 ; p++){
        updated_tab << mean_dph[et][p];
        if(p == size_mom_bin-2) updated_tab << "\n";
        else updated_tab << "\t";
      }
    }
    updated_tab << "\n";
    for(int et = 0 ; et < size_eta_bin-1 ; et++){
      for(int p = 0 ; p < size_mom_bin-1 ; p++){
        updated_tab << mean_des[et][p];
        if(p == size_mom_bin-2) updated_tab << "\n";
        else updated_tab << "\t";
      }
    }
    updated_tab.close();
  }
  // -------------------------------------------------------------
  // Manipulate Histogram Y-Axis
  better_yaxis(h1_dpp_v_p_et_bins,size_eta_bin-1);
  better_yaxis(h1_dth_v_p_et_bins,size_eta_bin-1);
  better_yaxis(h1_dph_v_p_et_bins,size_eta_bin-1);
  better_yaxis(h1_dpp_v_et_p_bins,size_mom_bin-1);
  better_yaxis(h1_dth_v_et_p_bins,size_mom_bin-1);
  better_yaxis(h1_dph_v_et_p_bins,size_mom_bin-1);

  // -------------------------------------------------------------
  // Plotting histograms
  TCanvas * c1 = new TCanvas("c1","c1",1300,900);
  c1 -> Divide(3,2);

  c1 -> cd(1); gPad -> SetRightMargin(0.04); gPad -> SetTopMargin(0.04); gPad ->SetLeftMargin(0.15); gPad -> SetLogy();
  h1_dpp_v_p_et_bins[0] -> Draw();
  for(int et = 0 ; et < size_eta_bin-1 ; et++) h1_dpp_v_p_et_bins[et] -> Draw("same");
  c1 -> cd(2); gPad -> SetRightMargin(0.04); gPad -> SetTopMargin(0.04); gPad ->SetLeftMargin(0.15); gPad -> SetLogy();
  h1_dth_v_p_et_bins[0] -> Draw();
  for(int et = 0 ; et < size_eta_bin-1 ; et++) h1_dth_v_p_et_bins[et] -> Draw("same");
  c1 -> cd(3); gPad -> SetRightMargin(0.04); gPad -> SetTopMargin(0.04); gPad ->SetLeftMargin(0.15); gPad -> SetLogy();
  h1_dph_v_p_et_bins[0] -> Draw();
  for(int et = 0 ; et < size_eta_bin-1 ; et++) h1_dph_v_p_et_bins[et] -> Draw("same");
  c1 -> cd(4); gPad -> SetRightMargin(0.04); gPad -> SetTopMargin(0.04); gPad ->SetLeftMargin(0.15); gPad -> SetLogy();
  h1_dpp_v_et_p_bins[0] -> Draw();
  for(int p = 0 ; p < size_mom_bin-1 ; p++) h1_dpp_v_et_p_bins[p] -> Draw("same");
  c1 -> cd(5); gPad -> SetRightMargin(0.04); gPad -> SetTopMargin(0.04); gPad ->SetLeftMargin(0.15); gPad -> SetLogy();
  h1_dth_v_et_p_bins[0] -> Draw();
  for(int p = 0 ; p < size_mom_bin-1 ; p++) h1_dth_v_et_p_bins[p] -> Draw("same");
  c1 -> cd(6); gPad -> SetRightMargin(0.04); gPad -> SetTopMargin(0.04); gPad ->SetLeftMargin(0.15); gPad -> SetLogy();
  h1_dph_v_et_p_bins[0] -> Draw();
  for(int p = 0 ; p < size_mom_bin-1 ; p++) h1_dph_v_et_p_bins[p] -> Draw("same");
  TLegend * leg1 = new TLegend(0.50,0.6,0.95,0.95);
  leg1 -> SetLineColor(0);
  leg1 -> SetFillStyle(0);
  for(int et = 0 ; et < size_eta_bin-1 ; et++) leg1 -> AddEntry(h1_dph_v_p_et_bins[et],Form("%.1f < |#eta| < %.1f",eta_bin[et],eta_bin[et+1]));

  if (N_Missing_Cut)
    c1 -> cd(2);
  else
    c1 -> cd(1);

  leg1 -> Draw("same");
  TLegend * leg2 = new TLegend(0.20,0.5,0.65,0.95);
  leg2 -> SetLineColor(0);
  leg2 -> SetFillStyle(0);
  for(int p = 0 ; p < size_mom_bin-1 ; p++) leg2 -> AddEntry(h1_dph_v_et_p_bins[p],Form("%.1f < P < %.1f GeV/c",mom_bin[p],mom_bin[p+1]));

  c1 -> cd(4);
  leg2 -> Draw("same");
  c1 -> Modified();
  c1 -> Update();

  // -------------------------------------------------------------
  // Saving fits to pdf
  for(int et = 0 ; et < size_eta_bin-1 ; et++){
    TString fname = out_pdf;
    if(et == 0) fname+="(";
    else if(et == size_eta_bin-2) fname+=")";
    c_fits_p [et] -> Print("fits/dpp_"+fname);
    c_fits_th[et] -> Print("fits/dth_"+fname);
    c_fits_ph[et] -> Print("fits/dph_"+fname);
  }

  // -------------------------------------------------------------
  // Saving histograms
  TFile * Fout = new TFile(outfile,"recreate");
  mom_response->Write();
  for(int et = 0 ; et < size_eta_bin-1 ; et++){
    h1_dpp_v_p_et_bins[et] -> Write(Form("h1_dpp_v_p_et_bins_%i",et));
    h1_dth_v_p_et_bins[et] -> Write(Form("h1_dth_v_p_et_bins_%i",et));
    h1_dph_v_p_et_bins[et] -> Write(Form("h1_dph_v_p_et_bins_%i",et));
  }
  for(int p = 0 ; p < size_mom_bin-1 ; p++){
    h1_dpp_v_et_p_bins[p] -> Write(Form("h1_dpp_v_et_p_bins_%i",p));
    h1_dth_v_et_p_bins[p] -> Write(Form("h1_dth_v_et_p_bins_%i",p));
    h1_dph_v_et_p_bins[p] -> Write(Form("h1_dph_v_et_p_bins_%i",p));
  }
  c1 -> Write("c1");
  TVT_eta_bin.Write("TVT_eta_bin");
  TVT_mom_bin.Write("TVT_mom_bin");
  Fout -> Close();

  // -------------------------------------------------------------
  // Saving plots
  c1 -> Print(out_pdf2);

  if(keep_plots)
    myapp -> Run();
  return 0;
}
// ============================================================================================================================================
void prettyTH1F( TH1F * h1 , int color , int marker , float min , float max ){
  int root_colors[8] = {52,56,66,71,78,90,94,100};
  h1 -> SetLineWidth(2);
  h1 -> SetLineColor(root_colors[color]);
  h1 -> SetMarkerStyle(marker);
  h1 -> SetMarkerColor(root_colors[color]);

  h1 -> GetXaxis() -> CenterTitle();
  h1 -> GetXaxis() -> SetNdivisions(107); // to draw less tick marks
  h1 -> GetYaxis() -> CenterTitle();
  h1 -> GetYaxis() -> SetNdivisions(107); // to draw less tick marks
}

void better_yaxis(TH1F ** h1_array,int array_size){

  //-----------------------------------------------
  //Determine the Y-axis min and max for TH1F Overlays
  float min = 999;
  float max = -999;

  for(int i = 0 ; i < array_size; i++){
    int min_bin = h1_array[i]->GetMinimumBin();
    int max_bin = h1_array[i]->GetMaximumBin();

    float temp_min = h1_array[i]->GetBinContent(min_bin)
      - h1_array[i]->GetBinError(min_bin);

    if (temp_min < 0)
      temp_min = h1_array[i]->GetBinContent(min_bin);

    float temp_max = h1_array[i]->GetBinContent(max_bin)
      + h1_array[i]->GetBinError(max_bin);

    if ((temp_min > 0) && (min > temp_min))
      min = temp_min;
    if (max < temp_max){
      max = temp_max;
    }
  }

  //-----------------------------------------------
  //Set the Y-axis Range
  for(int i = 0 ; i < array_size; i++){
    h1_array[i]->GetYaxis()->SetRangeUser(min,max);
  }
  //-----------------------------------------------
  //Labels
  h1_array[0]->GetYaxis()->SetTitleOffset(2.3);
  h1_array[0]->GetYaxis()->SetMoreLogLabels();
  return;
}

TF1 * double_gaus(TH1F *h1,float min1, float max1, float min2, float max2,TString type ,int et,int p){
  //First range is for narrow peak. Second range is full fit range desired.

  Double_t par[6];
  TF1 *G1 = new TF1 (Form("G1_%i_%i",et,p),"gaus",min1,max1);
  TF1 *G2 = new TF1 (Form("G2_%i_%i",et,p),"gaus",min2,max2);
  h1->Fit(G1,"R"); G1->GetParameters(&par[0]);
  h1->Fit(G2,"R"); G2->GetParameters(&par[3]);

  TF1 *double_gaus = new TF1(Form("f_gaus_%s_%i_%i",type.Data(),et,p),"gaus(0)+gaus(3)",min2,max2);
  double_gaus->SetParameters(par);
  return double_gaus;
}

float constituent_pT_threshold(float eta)
{
  // Minimum pT for B = 1.5 T (https://physdiv.jlab.org/DetectorMatrix/):

  // 100 MeV/c for -3.0 < eta < -2.5
  // 130 MeV/c for -2.5 < eta < -2.0
  // 70 MeV/c for -2.0 < eta < -1.5
  // 150 MeV/c for -1.5 < eta < -1.0
  //200 MeV/c for -1.0 < eta < 1.0
  // 150 MeV/c for 1.0 < eta < 1.5
  // 70 MeV/c for 1.5 < eta < 2.0
  // 130 MeV/c for 2.0 < eta < 2.5
  // 100 MeV/c for 2.5 < eta < 3.0

  eta = abs(eta);
  float eta_bins[6] = {3.5,2.5,2.0,1.5,1.0,0.0};
  float pT_threshold_array[5] = {0.1,0.13,0.70,0.15,0.2}; 
  float pT_threshold = 0;
  for (int i = 0; i < 5; i++)
    if ( (eta < eta_bins[i])&&(eta > eta_bins[i+1]) )
      pT_threshold = pT_threshold_array[i];

  return pT_threshold;
}

TGraphErrors ** Th1_to_TGraph(TH1F ** h1_array, int array_size)
{
  int root_colors[8] = {52,56,66,71,78,90,94,100};
  for(int i = 0 ; i < array_size; i++){
    TGraph *g = new TGraph();

    g->SetTitle(h1_array[i]->GetTitle());
    for (int i = 1; i <= h1_array[i]->GetNbinsX(); i++)
    {
      if (h1_array[i]->GetBinContent(i) == 0.0) continue;
      g->SetPoint(g->GetN(), h1_array[i]->GetBinCenter(i), h1_array[i]->GetBinContent(i));
      g->SetPoint(g->GetN(), h1_array[i]->GetBinCenter(i), h1_array[i]->GetBinError(i));
    }
    g->SetLineColor(root_colors[i]);
    g->SetMarkerStyle(8);
    g->SetMarkerColor(root_colors[i]);
  }
}
