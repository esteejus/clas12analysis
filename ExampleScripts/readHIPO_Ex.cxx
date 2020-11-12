#include <iostream>
#include "TH1.h"
#include "HipoChain.h"
#include "clas12reader.h"

using namespace clas12root;
//using namespace clas12root;
//
// "ROOT Script" entry point (the same name as the "filename's base").
//
// [bash/csh] root readHIPO_Ex.cxx
// [bash/csh] root readHIPO_Ex.cxx++
// root [0] .x readHIPO_Ex.cxx
// root [0] .x readHIPO_Ex.cxx++
//
void readHIPO_Ex(void)
{
  double emass = .511;
  HipoChain chain;
  chain.Add("/home/justin/analysis/JLabHIPO/cooked_TestHIPO/inc_011287.hipo");

  //loop over files
  auto c12 = chain.GetC12Reader();
  
  //create particles before loop
  TLorentzVector p4_elect1;
  TLorentzVector p4_elect2;

  //histograms
  TH1D *hmass = new TH1D("elec mass","Invariant mass electron",100,0,1);  

  c12->addExactPid(11,1);   //exactly 1 electrons
  c12->addExactPid(2212,1);   //exactly 1 proton
  //  c12->addAtLeastPid(11,2);   //exactly 2 electrons
  c12->addZeroOfRestPid(); //nothing else
  c12->useFTBased();        //and use the Pids from RECFT
  
  int counter = 0;

  //loop over all events in the file
  while(c12->next()==true)
    {
      c12 = chain.GetC12Reader();
      cout<<counter<<endl;
      counter++;
      
      if(c12->getDetParticles().empty())
	continue;

      auto elect1 = c12->getByID(11)[0];
      auto proton1 = c12->getByID(2212)[0];      

      p4_elect1.SetXYZM(elect1->par()->getPx(),elect1->par()->getPy(),elect1->par()->getPz(),emass);
      p4_elect2.SetXYZM(proton1->par()->getPx(),proton1->par()->getPy(),proton1->par()->getPz(),emass);
    

      //Fill histograms if electrons are in FD
      if(elect1->getRegion()==FD && proton1->getRegion()==FD)
	{
	  hmass->Fill(p4_elect1.M());

	}

    }

  hmass->Draw();




}

#if !defined(__CINT__) && !defined(__ACLIC__)
//
// "Standalone Application" entry point ("main").
//
//`root-config --cxx --cflags` -I/home/justin/clas12root/ -I/home/justin/clas12root/Clas12Banks/ -I/home/justin/clas12root/hipo4/ -o readHIPO_Ex readHIPO_Ex.cxx `root-config --libs`
// ./readHIPO_Ex
//
int main(int /*argc*/, char ** /*argv*/)
{
  readHIPO_Ex(); // just call the "ROOT Script"
  return 0;
}
#endif /* !defined(__CINT__) && !defined(__ACLIC__) */
