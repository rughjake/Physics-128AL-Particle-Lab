// The program to search for neutral pions.
// The neutral pion (pizero) decays into two photons
// The mass of the pizero is 0.135GeV
// Let us combine two photons to calculate the neutral pion mass.
//+
// File : pi0.cc
// Description : Analyze a data file
//
// Author : Ryosuke Itoh, IPNS, KEK
// Date : 1 - Mar - 2004
//-
// Modified by T.Nozaki  (2005/10/28)
// Modified by S.Nishida (2007/02/03)
// Modified by T.Nozaki (2009/07/01) Comments in English

// Modified by Jake Rugh (5/5/2020) Debugged, added histogram axis labels,
// added Lorentzian fit to histogram
#include <cstdio>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TClonesArray.h"
#include "TBrowser.h"
#include "TMath.h"

#include "BParticle.h"
#include "BEvent.h"
//enum SIMPLEPID {PHOTON, ELECTRON, PION, MUON, KAON, PROTON};

using namespace std;

// Quadratic background function
Double_t background(Double_t *x, Double_t *par) {
   return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

// Lorenzian Peak function
Double_t lorentzianPeak(Double_t *x, Double_t *par) {
  return (0.5*par[0]*par[1]/TMath::Pi()) /
    TMath::Max( 1.e-10,(x[0]-par[2])*(x[0]-par[2])
   + .25*par[1]*par[1]);
}

// Sum of background and peak function
Double_t fitFunction(Double_t *x, Double_t *par) {
  return background(x,par) + lorentzianPeak(x,&par[3]);
}

void analysis ( char* file, int maxevt = 0 ) {
 /* These are initialization codes.
     You can ignore them as a black box. */

  // Open a data file
  TFile *f = new TFile( file );

  // Obtain a pointer to a series of "event" data in the file
  TTree *t = (TTree*)f->Get("T");

  // Create a pointer to "BEvent" object where data are loaded from the file
  BEvent* event = new BEvent();

  // Obtain a branch for "BEvent" in the tree
  TBranch *branch = t->GetBranch("BEvent");

  // register created "BEvent" object to the branch to load event data
  branch->SetAddress(&event);

  
  /* Initialization codes finish here.*/


  // ---------------------------------------------------------------------- //

  /* We define the  histograms.*/

  // First, we define the name of the file where histograms will be saved. 
  // We define "histfile" as the histogram file name, as shown below. 
  TFile *hf = new TFile("histfile", "RECREATE" );

  // We define each histogram.  
  // The first one is the distribution of the number of particles in each event.
  // The longitudinal axis is divided into 50 bins between 0 and 50. 
  // The first argument in the blacket is the name of the histogram and 
  // the second one is its explanation.
  TH1F *h_nprt = new TH1F ( "h_nprt", "No. of particles", 50, 0, 50 );

  // The mass distribution obtained by combining two photons.
  // The longitudinal axis is divided into 100 bins between 0.0 and 1.0 GeV. 
  TH1F *h_pi0mass1 = new TH1F ( "h_pi0mass1", "pi0 mass with PID",
  				100, 0.0, 1.0 );
  h_pi0mass1->GetXaxis()->SetTitle("mass (GeV)");
  // The longitudinal axis is divided into 500 bins between 0.0 and 1.0 GeV. 
  TH1F *h_pi0mass2 = new TH1F ( "h_pi0mass2", "pi0 mass with PID",
				250, 0.0, 0.5 );
  h_pi0mass2->GetXaxis()->SetTitle("mass (GeV)");

  TH1F *h_pi0_LorFit = new TH1F ( "h_pi0_LorFit", "pi0 mass with fitted Lorentzian",
                                250, 0.0, 0.5 );
  h_pi0_LorFit->GetXaxis()->SetTitle("mass (GeV)");

  // ---------------------------------------------------------------------- //

  /* The analysis codes start here. */

  // Let us examine how many events exists in the data file.
  int nevent = (int)t->GetEntries();

  // Let us determine how many events should be analyzed.
  // If maxevent is not specified all the events in the data are analyzed.
  int nevent_process = nevent;
  if( maxevt>0 && maxevt<nevent ) nevent_process = maxevt;

  
  // Analysis is carried out event by event.
  // loop of events.
  for ( int i=0;i<nevent_process;i++ ) {
    //  Let us write an event number every 1000 events.
    // It helps for debugging when the program is dead.
    if ( i%1000 == 0 ) printf ( "***** Event %d\n", i );

    // An event is read.
    int nb = t->GetEntry ( i );

    // If the comment delimeter in the following sentence is removed, 
    // the detailed information of the first event is shown.
    // if ( i == 0 ) t->Show();

   /* -------- The main body of the analysis codes--------------- */

    // Obtain the number of particles in each event.
    int nparticles = event->NParticles();
    // It is stored in the histogram named h_nprt.
    h_nprt->Fill ( nparticles );

    // Particles aretreated as an array.
    TClonesArray& plist = *(event->GetParticleList());

    /* Let us combine two photons to calculate the neutral pion mass.
       When you want to change the program, you only need to modify the following part
       and the definition of the histogram. */

    // Loops to combine two particles.
    // When you combine two particles of the same type, the index k in the second for-loop
    // must begin with  k=j+1 to exclude double counting.
    for ( int j=0;j<nparticles;j++ ) {
      for ( int k=j+1;k<nparticles;k++ ) {
	// Selected two particles are named p1 and p2.
	BParticle* p1 = (BParticle*)plist[j];
	BParticle* p2 = (BParticle*)plist[k];

        // Charges of p1 and  p2 are checked.
        int charge1 = p1->charge();
        int charge2 = p2->charge();

	// Types of p1 and  p2 are checked.
	int pid1 = p1->pid();
	int pid2 = p2->pid();

        // Since the neutral pion decays into two photons
       // and the charge of the photon is zero,
       // check of the charge is actually not needed.       
        if( charge1==0 && charge2==0 ) {  // Both the first and second are neutral.
          if( pid1==PHOTON && pid2==PHOTON ) { // Both the first and second are photon.
            // Combine two particles to calculate the neutral pion mass.
            float px = p1->px()+p2->px(); // x component of momentum
            float py = p1->py()+p2->py(); // y component of momentum
            float pz = p1->pz()+p2->pz(); // z component of momentum
            float e = p1->e()+p2->e();    // Energy
 	    float mass = sqrt ( e*e-px*px-py*py-pz*pz ); // Calculate the combined mass.
	    h_pi0mass1->Fill ( mass );  // Store mass in h_pi0mass1.
	    h_pi0mass2->Fill ( mass );  // Store mass in h_pi0mass2.
            h_pi0_LorFit->Fill ( mass );
          }
        }
      }
    }
    //  Clear at the end of particle loop.
    event->Clear();
  }
  // Let us show how many events are processed.
  printf( "***** %d events processed. Exit.\n", nevent_process );


  double bin1 = h_pi0_LorFit->FindFirstBinAbove(h_pi0_LorFit->GetMaximum()/2);
  double bin2 = h_pi0_LorFit->FindLastBinAbove(h_pi0_LorFit->GetMaximum()/2);
  double gamma = h_pi0_LorFit->GetBinCenter(bin2) - h_pi0_LorFit->GetBinCenter(bin1);
  double x0 = h_pi0_LorFit->GetBinCenter(h_pi0_LorFit->GetMaximumBin());

  TF1 *L = new TF1("L",fitFunction,0.1,0.17,6);
  L->SetParameters(1,1,1,1,gamma,x0);
  h_pi0_LorFit->Fit("L");

//  h_pi0mass1->Draw("E");
//  h_pi0mass2->Draw("E");
  h_pi0_LorFit->Draw("E");
  gPad->Update();
  TPaveStats *st = (TPaveStats*)h_pi0_LorFit->FindObject("stats");
  gStyle->SetOptFit(0100);


  //  Let us store  all the histogram data in the file named histfile.
  hf->Write();
}

// The function showhist is used when the histograms are displayed.
void showhist() {
  TFile *f = new TFile ( "histfile" );
  TBrowser *t = new TBrowser;
}
