#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <filesystem>
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
