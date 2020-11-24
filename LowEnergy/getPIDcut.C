//=================================
//This is a quick macro to take slices in x-axis, fit a guassian to the y-distribution ...
//...and draw a TCutG cut contour around the range (mean - mult*sigma, mean + mult*sigma) for each x-bin
//to run change Input parameters and  type "root getPIDcut.C"


TGraph* DrawMean(TH1 *mean, TH1 *sigma)
{
  TGraph *g = new TGraph(mean->GetNbinsX());

  for(int i = 0; i < mean->GetNbinsX(); i++)
    g->SetPoint(i,mean->GetBinCenter(i),mean->GetBinContent(i));

  return g;
}

TCutG* DrawCut(TH1 *mean, TH1 *sigma, double mult)
{
  vector<double> x,y;
  for(int i = 0; i < mean->GetNbinsX(); i++)
    {
      if(mean->GetBinContent(i) < 1e-5)
	continue;

      x.push_back(mean->GetBinCenter(i));
      y.push_back(mean->GetBinContent(i) + sigma->GetBinContent(i)*mult);
    }


  int firstbin = 0;
  for(int i = 0; i < mean->GetNbinsX(); i++)
    {
      if(mean->GetBinContent(mean->GetNbinsX()-i-1) < 1e-5)
      continue;

      x.push_back(mean->GetBinCenter(mean->GetNbinsX()-i-1));
      y.push_back(mean->GetBinContent(mean->GetNbinsX()-i-1) - sigma->GetBinContent(sigma->GetNbinsX()-i-1)*mult);
    }

  //connect the dots
  x.push_back(x.front());
  y.push_back(y.front());

  
  return (new TCutG("cut",x.size(),x.data(),y.data()) );
}


void getPIDcut()
{
  //================
  //Input parameters

  TFile *inFile = TFile::Open("./hists/hist_merge.root");
  TString hist_name = "pid_p_ec";
  double mult = 2.; //multiple of sigma for cut
  vector<double> fit_range = {1.,4.}; //x-axis limits

  //===============

  TH2D *hist = (TH2D *)inFile->Get(hist_name);
  //root FitSlicesY makes histograms of all parameters for each x-bin, default is Guassian with norm, mean, sigma.
  hist->FitSlicesY(0,hist->GetXaxis()->FindBin(fit_range.at(0)),hist->GetXaxis()->FindBin(fit_range.at(1)));

  //Get these histograms
  TH1D *par1 = (TH1D *)inFile->Get((TString)hist->GetName()+"_1"); //par 1 is the mean
  TH1D *par2 = (TH1D *)inFile->Get((TString)hist->GetName()+"_2"); //par 2 is the sigma


  auto cut_g = DrawCut(par1,par2,mult); //Get the cut ranging from (mean - sigma * mult , mean + sigma * mult)
  auto mean_g = DrawMean(par1,par2); //Get the mean value vs x-axis variable

  TCanvas *cvs = new TCanvas(1);
  mean_g->SetLineColor(2);
  hist->Draw("colz");
  mean_g->Draw("same");
  cut_g->Draw("same");
  
}
