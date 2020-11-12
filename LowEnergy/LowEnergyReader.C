#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TCutG.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"

using namespace clas12;


void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
	     rp->par()->getPz(),p4.M());

}

void LowEnergyReader(){
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();


  /////////////////////////////////////
  //ignore this just getting file name!
  TString inputFile;
  TString outputFile;

  for(Int_t i=1;i<gApplication->Argc();i++){
    TString opt=gApplication->Argv(i);
    if((opt.Contains(".hipo"))){
      inputFile=opt(5,opt.Sizeof());
    }
  }
  if(inputFile==TString())  {
    std::cout << " *** please provide a file name..." << std::endl;
    exit(0);
  }
  /////////////////////////////////////

   
  cout<<"Analysing hipo file "<<inputFile<<endl;

  TChain fake("hipo");
  fake.Add(inputFile.Data());
  //get the hipo data
  //   reader.open(inputFile.Data());
  auto files=fake.GetListOfFiles();

  //some particles
  //  clas12databases::SetCCDBLocalConnection("/home/justin/ccdb/ccdb.sqlite");
  //  clas12databases::SetQADBConnection("/home/justin/clasqaDB/qaDB.json");
  //  clas12databases::SetRCDBRootConnection("/home/justin/rcdb/rcdb.root");

  auto pdg_db=TDatabasePDG::Instance();
  TLorentzVector beam(0,0,4.17179,4.17179);
  TLorentzVector target(0,0,0, 1.8917 + .013 );
  TLorentzVector el(0,0,0,pdg_db->GetParticle(11)->Mass());
  TLorentzVector pr(0,0,0,pdg_db->GetParticle(2212)->Mass());

  auto *pid_e_ec = new TH2D("pid_p_ec","EC/p vs e electron",100,0,5,100,0,.6);
  auto *pid_p_fd = new TH2D("pid_p_fd","FD TOF vs p PID proton",500,0,4,500,0,1.3);
  auto *pid_p_cd = new TH2D("pid_p_cd","CD TOF vs p PID proton",200,0,4,200,0,1.3);

  auto *chi2_p = new TH1D("chi2_p","Chi2 proton",100,-20,20);
  auto *chi2_e = new TH1D("chi2_e","Chi2 electron",100,-20,20);

  auto *hmiss=new TH1F("missM","missM",200,-1,4);
  auto *missm_pmiss = new TH2F("missm_pmiss","Missing mass vs. pmiss",200,-1,4,100,0,4);
  auto *missm_xb = new TH2F("missm_xb","Missing mass vs. x_{B}",200,-1,4,100,0,4);
  auto *theta_pq_ratio = new TH2D("theta_pq_ratio","Theata pq vs p/q",50,0,1,50,0,50);

  auto *xbor=new TH1F("xBor","xBor",200,-1,10);   

  gBenchmark->Start("timer");
  int counter=0;
 
  TFile *cutf = TFile::Open("cuts.root");
  TCutG *elec_cut = (TCutG *)cutf->Get("elec");
  TCutG *ftof_p = (TCutG *)cutf->Get("ftof_p");
  TCutG *ctof_p = (TCutG *)cutf->Get("ctof_p");
  
  for(Int_t i=0;i<files->GetEntries();i++){
    //create the event reader
    clas12reader c12(files->At(i)->GetTitle(),{0});
    //    clas12databases db;
    //    c12.connectDataBases(&db);
    
    //  clas12reader c12(files->At(i)->GetTitle(),{0});//add tags {tag1,tag2,tag3,...}
      
    //Add some event Pid based selections
    //////////c12.AddAtLeastPid(211,1); //at least 1 pi+
    //c12.addExactPid(11,1);    //exactly 1 electron
    //c12.addExactPid(211,1);    //exactly 1 pi+
    //c12.addExactPid(-211,1);    //exactly 1 pi-
    //c12.addExactPid(2212,1);    //exactly 1 proton
    //c12.addExactPid(22,2);    //exactly 2 gamma
    //////c12.addZeroOfRestPid();  //nothing else
    //////c12.useFTBased(); //and use the Pids from RECFT

    //can also access the integrated current at this point
    //c12.scalerReader();//must call this first
    //c12.getRunBeamCharge();
    //    c12.getRcdbVals();

    //c12.setEntries(1E3); //only process 1E3 events per file
    while(c12.next()==true){
      //can get an estimate of the beam current to this event
      //c12.getCurrApproxCharge();//if called c12.scalerReader();
	
      //c12.event()->getStartTime();

	
      //Loop over all particles to see how to access detector info.
    
      for(auto& p : c12.getDetParticles()){
	//  get predefined selected information
	p->getTime();
	p->getDetEnergy();
	p->getDeltaEnergy();

	//check trigger bits
	//	 if(c12.checkTriggerBit(25)) cout<<"MesonExTrigger"<<endl;
	//	 else cout<<"NOT"<<endl;

	// get any detector information (if exists for this particle)
	// there should be a get function for any entry in the bank
	switch(p->getRegion()) {
	case FD :
	  p->cal(PCAL)->getEnergy();
	  p->cal(ECIN)->getEnergy();
	  p->cal(ECOUT)->getEnergy();
	  p->sci(FTOF1A)->getEnergy();
	  p->sci(FTOF1B)->getEnergy();
	  p->sci(FTOF2)->getEnergy();
	  p->trk(DC)->getSector();
	  p->trk(DC)->getChi2();
	  p->che(HTCC)->getNphe();
	  p->che(LTCC)->getNphe();
	  //trajectories
	  p->traj(LTCC)->getX();
	  // p->traj(DC,DC1)->getCx();; //First layer of DC, hipo4
	  break;
	case FT :
	  p->ft(FTCAL)->getEnergy();
	  p->ft(FTHODO)->getEnergy();
	  break;
	case CD:
	  p->sci(CTOF)->getEnergy();
	  p->sci(CND)->getEnergy();
	  break;
	}
	//   covariance matrix (comment in to see!)
	// p->covmat()->print();
	p->cmat();
      }

      // get particles by type
      auto electrons=c12.getByID(11);
      auto protons=c12.getByID(2212);
      //       auto gammas=c12.getByID(22);
      //       auto pips=c12.getByID(211);
      //       auto pims=c12.getByID(-211);
       
      //      if(protons.size()==1)
	//	cout<< protons[0]->getRegion()<<endl;

      if(electrons.size()==1 && protons.size()==1 ){
       

	//       	if(protons[0]->getRegion()==FD && !ftof_p ->IsInside(protons[0]->par()->getP(),protons[0]->par()->getBeta()) )
	//	  continue;

	//	if( protons[0]->getRegion()==CD && ctof_p ->IsInside(protons[0]->par()->getP(),protons[0]->par()->getBeta() )
	   //	    continue;
	    
	//	if( electronss[0]->getRegion()==FD && elec_cut ->IsInside(electrons[0]->par()->getP(),electrons[0]->cal(PCAL)->getEnergy()/electrons[0]-par()->getP() )
	   //	    continue;

	if(protons[0]->getRegion()==FD)
	  pid_p_fd -> Fill(protons[0]->par()->getP(),protons[0]->par()->getBeta() );

	if(protons[0]->getRegion()==CD)
	  pid_p_cd -> Fill(protons[0]->par()->getP(),protons[0]->par()->getBeta() );

	if(electrons[0]->getRegion()==FD)
	  pid_e_ec -> Fill(electrons[0]->par()->getP(), electrons[0]->cal(PCAL)->getEnergy()/electrons[0]->par()->getP() );
	  
	//	chi2_p->Fill(protons[0]->par()->getChi2Pid() );
	//	chi2_e->Fill(electrons[0]->par()->getChi2Pid() );
	
	// set the particle momentum
	SetLorentzVector(el,electrons[0]);
	SetLorentzVector(pr,protons[0]);
	
	TLorentzVector miss=beam+target-el-pr; //missing 4-vector
	TLorentzVector q = beam - el;          //photon  4-vector
	double x_b = -q.M2()/(2*pdg_db->GetParticle(2212)->Mass()*q.E() ); //x-borken
	double theta = pr.Vect().Angle(q.Vect()) * TMath::RadToDeg();  //angle between vectors p_miss and q


	bool pcut_fd = (protons[0]->getRegion()==FD && ftof_p ->IsInside(protons[0]->par()->getP(),protons[0]->par()->getBeta()) );
	bool pcut_cd = (protons[0]->getRegion()==CD && ctof_p ->IsInside(protons[0]->par()->getP(),protons[0]->par()->getBeta() ) );
	bool ecut = (electrons[0]->getRegion()==FD && elec_cut ->IsInside(electrons[0]->par()->getP(),electrons[0]->cal(PCAL)->getEnergy()/electrons[0]->par()->getP()) );


	if( (pcut_fd || pcut_cd) && ecut)
	  {
	    hmiss->Fill(miss.M());
	    missm_pmiss->Fill(miss.M(),miss.P());
	    missm_xb->Fill(miss.M(),x_b);
	    theta_pq_ratio->Fill(pr.Vect().Mag()/q.Vect().Mag(),theta);
	  }

	//could also get particle time etc. here too
	//Double_t eTime=electrons[0]->sci(FTOF1A)->getTime();
      }
    
       
      counter++;
    }

  }
  gBenchmark->Stop("timer");
  gBenchmark->Print("timer");

  theta_pq_ratio->GetYaxis()->SetTitle("#theta_{pq}");
  theta_pq_ratio->GetXaxis()->SetTitle("|p|/|q|");
  missm_pmiss->GetYaxis()->SetTitle("Missing mass (GeV)");
  missm_pmiss->GetXaxis()->SetTitle("P_{miss} (GeV)");
  pid_e_ec->GetYaxis()->SetTitle("E_{loss}/p");
  pid_e_ec->GetXaxis()->SetTitle("p (GeV)");				
  pid_p_fd->GetYaxis()->SetTitle("FD TOF");
  pid_p_fd->GetXaxis()->SetTitle("p_{proton} (GeV)");				
  pid_p_cd->GetYaxis()->SetTitle("CD TOF");
  pid_p_cd->GetXaxis()->SetTitle("p_{proton} (GeV)");				

  theta_pq_ratio->GetYaxis()->CenterTitle();
  theta_pq_ratio->GetXaxis()->CenterTitle();
  missm_pmiss->GetYaxis()->CenterTitle();
  missm_pmiss->GetXaxis()->CenterTitle();
  pid_e_ec->GetYaxis()->CenterTitle();
  pid_e_ec->GetXaxis()->CenterTitle();
  pid_p_fd->GetYaxis()->CenterTitle();
  pid_p_fd->GetXaxis()->CenterTitle();
  pid_p_cd->GetYaxis()->CenterTitle();
  pid_p_cd->GetXaxis()->CenterTitle();

  TCanvas* can = new TCanvas();
  can->Divide(2,2);
  can->cd(1);
  missm_pmiss->Draw("colz");
  can->cd(2);
  missm_xb->Draw("colz");
  can->cd(3);
  theta_pq_ratio->Draw("colz");
  
  TCanvas* can1 = new TCanvas();
  pid_e_ec->Draw("colz");
  /*
  can1->Divide(2,2);
  can1->cd(1);
  pid_p_fd->Draw("colz");
  can1->cd(2);
  pid_p_cd->Draw("colz");
  can1->cd(3);
  pid_e_ec->Draw("colz");
  */
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<counter<< " s\n";

}
