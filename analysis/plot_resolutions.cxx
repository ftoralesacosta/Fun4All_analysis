#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <filesystem>

// Root includes
#include "TRint.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TVectorT.h"
#include "TGraphErrors.h"

namespace fs = std::filesystem;
using namespace std;

// Forward-declaring functions
void prettyTH1F( TH1F * h1 , int color , int marker , float min , float max );
int idx_from_vector( double value , TVectorT<double> * vec );
void prettyTH1( TH1F * h1 , int color , int marker , float min , float max );
double sq(double x){ return x*x;}
// ============================================================================================================================================
int main(int argc, char ** argv) {

#ifdef WITHRINT
	TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
#else
	TApplication *myapp = new TApplication("myapp",0,0);
#endif

	//gStyle->SetErrorX(0.0001);
	gStyle->SetTitleSize(0.08,"t");
	// ------------------------------------------------------------------------------
	// List paths to files that will be loaded
	TString fnames[] = {
	  "../output/output_mom_res_combined_electrons+jetssigma_eta_5_p_5_.root"
	};
	// #######################################################################################################################################
	// YOU SHOULDN'T NEED TO MODIFY ANYTHING IN THE BLOCK OF CODE BELOW AND UNTIL AFTER THE NEXT LINE WITH ###...
	const int size_loaded = sizeof(fnames)/sizeof(*fnames);
	// ------------------------------------------------------------------------------
	// Preparing variables that will later on be filled from root files
	TVectorT<double> ** TVT_eta_bin = new TVectorT<double>*[size_loaded];
	TVectorT<double> ** TVT_mom_bin = new TVectorT<double>*[size_loaded];
	int num_eta_bin[size_loaded] = {0};
	int num_mom_bin[size_loaded] = {0};
	// ------------------------------------------------------------------------------
	// Loading root files and info therein
	TFile ** Fin = new TFile*[size_loaded];

	TH1F *** h1_dpp_v_p_et_bins = new TH1F ** [size_loaded];
	TH1F *** h1_dth_v_p_et_bins = new TH1F ** [size_loaded];
	TH1F *** h1_dph_v_p_et_bins = new TH1F ** [size_loaded];
	TH1F *** h1_dpp_v_et_p_bins = new TH1F ** [size_loaded];
	TH1F *** h1_dth_v_et_p_bins = new TH1F ** [size_loaded];
	TH1F *** h1_dph_v_et_p_bins = new TH1F ** [size_loaded];

	for(int f = 0 ; f < size_loaded ; f++){
		ifstream fin;
		fin.open(fnames[f]);
		if(!fin){ cout << "\033[1;31mCouldn't find input file '" << fnames[f] << "'. Bailing out!\033[0m" << endl; exit(0);}
		fin.close();
		cout<<__LINE__<<endl;
		Fin[f] = new TFile(fnames[f]);

		TVT_eta_bin[f] = (TVectorT<double> *) Fin[f] -> Get("TVT_eta_bin");
		TVT_mom_bin[f] = (TVectorT<double> *) Fin[f] -> Get("TVT_mom_bin");
		cout<<__LINE__<<endl;
		num_eta_bin[f] = (*TVT_eta_bin[f]).GetNoElements()-1;
		num_mom_bin[f] = (*TVT_mom_bin[f]).GetNoElements()-1;

                h1_dpp_v_p_et_bins[f] = new TH1F * [num_eta_bin[f]];
                h1_dth_v_p_et_bins[f] = new TH1F * [num_eta_bin[f]];
                h1_dph_v_p_et_bins[f] = new TH1F * [num_eta_bin[f]];
                h1_dpp_v_et_p_bins[f] = new TH1F * [num_mom_bin[f]];
                h1_dth_v_et_p_bins[f] = new TH1F * [num_mom_bin[f]];
                h1_dph_v_et_p_bins[f] = new TH1F * [num_mom_bin[f]];
		
                for(int et = 0 ; et < num_eta_bin[f] ; et++){
		  cout<<__LINE__<<endl;
		  h1_dpp_v_p_et_bins[f][et] = (TH1F*) Fin[f] -> Get(Form("h1_dpp_v_p_et_bins_%i",et));
                        h1_dth_v_p_et_bins[f][et] = (TH1F*) Fin[f] -> Get(Form("h1_dth_v_p_et_bins_%i",et));
                        h1_dph_v_p_et_bins[f][et] = (TH1F*) Fin[f] -> Get(Form("h1_dph_v_p_et_bins_%i",et));
                }
		
                for(int p = 0 ; p < num_mom_bin[f] ; p++){
		  cout<<__LINE__<<endl;
                        h1_dpp_v_et_p_bins[f][p ] = (TH1F*) Fin[f] -> Get(Form("h1_dpp_v_et_p_bins_%i",p ));
                        h1_dth_v_et_p_bins[f][p ] = (TH1F*) Fin[f] -> Get(Form("h1_dth_v_et_p_bins_%i",p ));
                        h1_dph_v_et_p_bins[f][p ] = (TH1F*) Fin[f] -> Get(Form("h1_dph_v_et_p_bins_%i",p ));
                }
        }
		
	cout<<__LINE__<<endl;	

	cout << "\neta bin boundaries:\n"; for(int et = 0 ; et < num_eta_bin[0]+1 ; et++) cout << (*TVT_eta_bin[0])[et] << ", "; cout << "\n";
	cout << "\np bin boundaries:\n"  ; for(int p  = 0 ; p  < num_mom_bin[0]+1 ; p ++) cout << (*TVT_mom_bin[0])[ p] << ", "; cout << "\n\n";

	// #######################################################################################################################################
        // EDIT THE CODE BELOW DEPENDING ON WHAT YOU WANT TO PLOT
	// ------------------------------------------------------------------------------
	double max[] = {29,13,7,3,2,2,2,2,2,2,2,2,3,7,13,29};
	// Editing histograms

	for(int f = 0 ; f < size_loaded ; f++){
	  for(int i = 0 ; i < num_eta_bin[0] ; i++){ 
	        h1_dpp_v_p_et_bins[f][i] -> SetTitle(Form("%.1f < #eta < %.1f",(*TVT_eta_bin[0])[i],(*TVT_eta_bin[0])[i+1]));
                h1_dpp_v_p_et_bins[f][i] -> SetTitle(Form("%.1f < #eta < %.1f",(*TVT_eta_bin[0])[i],(*TVT_eta_bin[0])[i+1]));
                h1_dpp_v_p_et_bins[f][i] -> SetTitle(Form("%.1f < #eta < %.1f",(*TVT_eta_bin[0])[i],(*TVT_eta_bin[0])[i+1]));
                h1_dpp_v_p_et_bins[f][i] -> SetTitle(Form("%.1f < #eta < %.1f",(*TVT_eta_bin[0])[i],(*TVT_eta_bin[0])[i+1])); 

		prettyTH1( h1_dpp_v_p_et_bins[f][i] ,  2 , 24 , 999 , max[i] );
		prettyTH1( h1_dpp_v_p_et_bins[f][i] , 96 , 20 , 999 , max[i] );
		prettyTH1( h1_dpp_v_p_et_bins[f][i] ,  4 , 24 , 999 , max[i] );
		prettyTH1( h1_dpp_v_p_et_bins[f][i] , 62 , 20 , 999 , max[i] );
        }

        TCanvas * c1 = new TCanvas("c1","c1",1300,900);
        c1 -> Divide(3,2);

        c1 -> cd(1); gPad -> SetRightMargin(0.04); gPad -> SetTopMargin(0.04); gPad ->SetLeftMargin(0.15); gPad -> SetLogy();
        h1_dpp_v_p_et_bins[f][0] -> Draw();
        for(int et = 0 ; et < num_eta_bin[f] ; et++) h1_dpp_v_p_et_bins[f][et] -> Draw("same");
        c1 -> cd(2); gPad -> SetRightMargin(0.04); gPad -> SetTopMargin(0.04); gPad ->SetLeftMargin(0.15); gPad -> SetLogy();
        h1_dth_v_p_et_bins[f][0] -> Draw();
        for(int et = 0 ; et < num_eta_bin[f] ; et++) h1_dth_v_p_et_bins[f][et] -> Draw("same");
        c1 -> cd(3); gPad -> SetRightMargin(0.04); gPad -> SetTopMargin(0.04); gPad ->SetLeftMargin(0.15); gPad -> SetLogy();
        h1_dph_v_p_et_bins[f][0] -> Draw();
        for(int et = 0 ; et < num_eta_bin[f] ; et++) h1_dph_v_p_et_bins[f][et] -> Draw("same");
        c1 -> cd(4); gPad -> SetRightMargin(0.04); gPad -> SetTopMargin(0.04); gPad ->SetLeftMargin(0.15); gPad -> SetLogy();
        h1_dpp_v_et_p_bins[f][0] -> Draw();
        for(int p = 0 ; p < num_mom_bin[f] ; p++) h1_dpp_v_et_p_bins[f][p] -> Draw("same");
        c1 -> cd(5); gPad -> SetRightMargin(0.04); gPad -> SetTopMargin(0.04); gPad ->SetLeftMargin(0.15); gPad -> SetLogy();
        h1_dth_v_et_p_bins[f][0] -> Draw();
        for(int p = 0 ; p < num_mom_bin[f] ; p++) h1_dth_v_et_p_bins[f][p] -> Draw("same");
        c1 -> cd(6); gPad -> SetRightMargin(0.04); gPad -> SetTopMargin(0.04); gPad ->SetLeftMargin(0.15); gPad -> SetLogy();
        h1_dph_v_et_p_bins[f][0] -> Draw();
        for(int p = 0 ; p < num_mom_bin[f] ; p++) h1_dph_v_et_p_bins[f][p] -> Draw("same");
        TLegend * leg1 = new TLegend(0.50,0.6,0.95,0.95);
        leg1 -> SetLineColor(0);
        leg1 -> SetFillStyle(0);
        for(int et = 0 ; et < num_eta_bin[f] ; et++) leg1 -> AddEntry(h1_dph_v_p_et_bins[f][et],Form("%.1f < |#eta| < %.1f",(*TVT_eta_bin[0])[et],(*TVT_eta_bin[0])[et+1]));

	//        if (N_Missing_Cut)
          c1 -> cd(2);
	  // else
          // c1 -> cd(1);

        leg1 -> Draw("same");
        TLegend * leg2 = new TLegend(0.20,0.5,0.65,0.95);
        leg2 -> SetLineColor(0);
        leg2 -> SetFillStyle(0);
        for(int p = 0 ; p < num_mom_bin[f] ; p++) leg2 -> AddEntry(h1_dph_v_et_p_bins[f][p],Form("%.1f < P < %.1f GeV/c",(*TVT_mom_bin[0])[p],(*TVT_mom_bin[0])[p]));

        c1 -> cd(4);
        leg2 -> Draw("same");
        c1 -> Modified();
        c1 -> Update();
        myapp -> Run();
	return 0;
 }
}


// ============================================================================================================================================
void prettyTH1F( TH1F * h1 , int color , int marker , float min , float max ){
	h1 -> SetLineWidth(2);
	h1 -> SetLineColor(color);
	h1 -> SetMarkerStyle(marker);
	h1 -> SetMarkerColor(color);

	if(min!=999) h1 -> SetMinimum(min);
	if(max!=999) h1 -> SetMaximum(max);

	h1 -> GetXaxis() -> CenterTitle();
	h1 -> GetXaxis() -> SetNdivisions(108); // to draw less tick marks
	h1 -> GetYaxis() -> CenterTitle();
	h1 -> GetYaxis() -> SetNdivisions(108); // to draw less tick marks

	h1 -> SetMinimum(0.001);
}
// ============================================================================================================================================
int idx_from_vector( double value , TVectorT<double> * vec ){
	int size_vec = (*vec).GetNoElements();
	double diff = 999.;
	int idx = -1;
	for(int i = 0 ; i < size_vec ; i++){
		if(abs((*vec)[i] - value)<diff){
			diff = abs((*vec)[i] - value);
			idx = i;
		}
	}
	return idx;
}
// ============================================================================================================================================
void prettyTH1( TH1F * h1 , int color , int marker , float min , float max ){
	h1->SetLineColor(color);
	h1->SetMarkerColor(color);
	h1->SetMarkerStyle(marker);
	h1->SetMarkerSize(1.4);
	h1->SetLineWidth(1);

	h1->GetXaxis()->SetNdivisions(108);
	h1->GetXaxis()->SetTitleSize(0.06);
	h1->GetXaxis()->SetLabelSize(0.06);
	h1->GetXaxis()->CenterTitle();

	h1->GetYaxis()->SetNdivisions(108);
	h1->GetYaxis()->SetTitleSize(0.06);
	h1->GetYaxis()->SetLabelSize(0.06);
	h1->GetYaxis()->CenterTitle();

	if(min!=999) h1 -> SetMinimum(min);
        if(max!=999) h1 -> SetMaximum(max);
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
