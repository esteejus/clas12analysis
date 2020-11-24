

void plotResults(TObjArray par, TF1 *f1, TH2D *h, vector<double> fit_range)
{

  double mult = 1; //how many sigma to integrate for signal/bkg ratio
  int st_bin  = h->GetYaxis()->FindBin(fit_range.at(0));
  int end_bin = h->GetYaxis()->FindBin(fit_range.at(1));
  int num = end_bin-st_bin;
  int dim = round(sqrt(num));

  TH1D *slice[num];
  TCanvas *c1 = new TCanvas("c1","c1",2000,2000);
  c1->DivideSquare(num,0,0);

  TF1 sig("sig","[0]*TMath::Gaus(x,[1],[2])",-.3,.3);
  TF1 bkg("sig","[0]*TMath::Gaus(x,[1],[2])",-.3,.3);
  TF1 gaus2("gaus2","[0]*TMath::Gaus(x,[1],[2]) + [3]*TMath::Gaus(x,[4],[5])",-.3,.3);
  gaus2.SetParameters(1.5e3,0,.05,1e3,2e-1,.3);

  sig.SetLineColor(2);
  bkg.SetLineColor(4);
  gaus2.SetLineColor(1);

  vector<double> pmiss_slice, sig_bkg_r, miss_mean, miss_resol;
  
  for(int i = 0; i < num; i++)
    {
      //      sig.SetParameters( ((TH1 *)par[0])->GetBinContent(st_bin + i), ((TH1 *)par[1])->GetBinContent(st_bin + i), ((TH1 *)par[2])->GetBinContent(st_bin + i));

      cout<<((TH1 *)par[0])->GetBinContent(st_bin + i)<<endl;
      slice[i] = h->ProjectionX(Form("slice_%d",i),st_bin + i,st_bin+1+i);
      c1->cd(i+1);
      c1->cd(i+1)->SetLogy();
      slice[i]->Fit(&gaus2,"","",-.3,.3);
      sig.SetParameters(gaus2.GetParameter(0),gaus2.GetParameter(1),gaus2.GetParameter(2));
      bkg.SetParameters(gaus2.GetParameter(3),gaus2.GetParameter(4),gaus2.GetParameter(5));
      cout<<"sig norm "<<gaus2.GetParameter(0)<<endl;
      slice[i]->GetXaxis()->SetRangeUser(-.5,.5);
      //      slice[i]->GetYaxis()->SetRangeUser(.1,2e3);
      //      slice[i]->GetYaxis()->SetRangeUser(10,2e4);

      //      slice[i]->SetTitle("");
      slice[i]->GetYaxis()->SetTitle("Counts");
      slice[i]->GetYaxis()->CenterTitle();
      slice[i]->GetXaxis()->SetTitle("Missing mass (GeV)");
      slice[i]->GetXaxis()->CenterTitle();


     
      slice[i]->GetFunction("gaus2")->SetBit(TF1::kNotDraw);
      slice[i]->Draw();
      gaus2.DrawCopy("same");
      sig.DrawCopy("same");
      bkg.DrawCopy("same");

      pmiss_slice.push_back(h->GetYaxis()->GetBinCenter(st_bin + i));
      cout<<"bin center "<<h->GetYaxis()->GetBinCenter(st_bin + i)<<endl;
      miss_mean.push_back(gaus2.GetParameter(1));
      miss_resol.push_back(abs(gaus2.GetParameter(2))*1000);
      sig_bkg_r.push_back(sig.Integral(gaus2.GetParameter(1)-mult*gaus2.GetParameter(2),gaus2.GetParameter(1)+mult*gaus2.GetParameter(2)) / bkg.Integral(gaus2.GetParameter(1)-mult*gaus2.GetParameter(2),gaus2.GetParameter(1)+mult*gaus2.GetParameter(2)) );
    }

  TGraph *g_miss_mean = new TGraph(pmiss_slice.size(),pmiss_slice.data(),miss_mean.data());
  TGraph *g_miss_resol = new TGraph(pmiss_slice.size(),pmiss_slice.data(),miss_resol.data());
  TGraph *g_sig_bkg = new TGraph(pmiss_slice.size(),pmiss_slice.data(),sig_bkg_r.data());
  g_miss_mean->SetTitle("Missing Mass Mean");
  g_miss_resol->SetTitle("Missing Mass Resolution (1#sigma)");
  g_sig_bkg->SetTitle("Signal/Background  (1#sigma around Sig.)");

  g_miss_mean->GetYaxis()->SetTitle("Missing mass mean - Neutron mass (GeV)");
  g_miss_mean->GetXaxis()->SetTitle("Missing momentum (GeV)");
  g_miss_resol->GetYaxis()->SetTitle("Missing Mass resolution (MeV)");
  g_miss_resol->GetXaxis()->SetTitle("Missing momentum (GeV)");
  g_sig_bkg->GetYaxis()->SetTitle("Signal to Bkg ratio");
  g_sig_bkg->GetXaxis()->SetTitle("Missing momentum (GeV)");

  g_miss_mean->GetYaxis()->CenterTitle();
  g_miss_mean->GetXaxis()->CenterTitle();
  g_miss_resol->GetYaxis()->CenterTitle();
  g_miss_resol->GetXaxis()->CenterTitle();
  g_sig_bkg->GetYaxis()->CenterTitle();
  g_sig_bkg->GetXaxis()->CenterTitle();

  g_miss_mean->SetMarkerStyle(20);
  g_miss_resol->SetMarkerStyle(20);
  g_sig_bkg->SetMarkerStyle(20);

  
  TCanvas *c2 = new TCanvas("c2","c2",1000,1000);
  c2->Divide(3,0,.02,.02);
  c2->cd(1);
  g_miss_mean->Draw("APO");
  c2->cd(2);
  g_miss_resol->Draw("APO");
  c2->cd(3);
  g_sig_bkg->Draw("APO");

  c1->SaveAs("fits.png");
  c2->SaveAs("resolution.png");
  return;
}



void fitMissMass()
{
  gStyle->SetOptStat(0);
  
  TFile *f = TFile::Open("./hist_merge_xb.root");
  TH2D *h = (TH2D *)f->Get("missm_pmiss");
  h->RebinY();
  //  h->RebinY();
  //  h->RebinY();
  
  vector<double> fit_range = {.4,2.};
  TF1 *gaus2 = new TF1("gaus2","[0]*TMath::Gaus(x,[1],[2]) + [3]*TMath::Gaus(x,[4],[5])",.6,1.2);
  gaus2->SetParameters(2.8e3,.94,.09,5.1e5,3.11,.6);

  TObjArray aSlices;
  h->FitSlicesX(gaus2,h->GetYaxis()->FindBin(fit_range.at(0)),h->GetYaxis()->FindBin(fit_range.at(1)),0,"Q",&aSlices);
  plotResults(aSlices,gaus2,h,fit_range);


}
