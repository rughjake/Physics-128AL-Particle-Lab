// File : D_meson.cc
// Author : Jake Rugh (5/5/2020)
//
// Edited from : phi_meson.cc

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

  // The mass distribution obtained by combining pi-plus and pi-minus.
  // The longitudinal axis is divided into 400 bins between 0.0 and 2.0 GeV. 
  TH1F *h_Dmass = new TH1F ( "h_Dmass", "D meson mass with PID", 400, 1.7, 2.0 );
  h_Dmass->GetXaxis()->SetTitle("mass (GeV)");

  TH1F *h_D_LorFit = new TH1F ( "h_D_LorFit", "D meson mass with PID", 400, 1.7, 2.0 );
  h_D_LorFit->GetXaxis()->SetTitle("mass (GeV)");

  // ---------------------------------------------------------------------- //

   /* The analysis codes start here. */

  // Let us examine how many events exists in the data file.
  int nevent = (int)t->GetEntries();

  // Let us determine how many events should be analyzed.
  // If maxevent is not specified all the events in the data are analyzed.
  int nevent_process = nevent;
  if( maxevt>0 && maxevt<nevent ) nevent_process = maxevt;

  /// Analysis is carried out event by event.
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

    /*-------- The main body of the analysis codes--------------- */

    // Obtain the number of particles in each event.
    int nparticles = event->NParticles();
    // It is stored in the histogram named h_nprt.
    h_nprt->Fill ( nparticles );

    // Particles aretreated as an array.
    TClonesArray& plist = *(event->GetParticleList());

    /* Let us combine pi-plus and pi-minus to calculate the neutral kaon mass.
       When you want to change the program, you only need to modify the following part
       and the definition of the histogram. */


    // Loops to combine two particles.
    // We do not adopt a code which makes a point of efficiency, but
    // adopt a code which makes a point of easiness in understanding. 
    // In the code which makes point of efficiency, k begins with k=j+1,
    // but in the present code k begins with k=0.
    // In the former code requirements on charges are 
    // if(charge1==-charge2 && charge1!=0)
    
    for ( int j=0;j<nparticles;j++ ) {
      for ( int k=0;k<nparticles;k++ ) {
        // Selected two particles are named p1 and p2.
        BParticle* p1 = (BParticle*)plist[j];
        BParticle* p2 = (BParticle*)plist[k];

        // Charges of p1 and  p2 are checked.
        int charge1 = p1->charge();
        int charge2 = p2->charge();

       // Types of p1 and  p2 are checked.
        int pid1 = p1->pid();
        int pid2 = p2->pid();

        // Let us combine pi-plus and pi-minus.
        if( charge1==1 && charge2==-1  ) {  // The first is positive, the second is negative.
          if( (pid1==PION && pid2==KAON) || (pid1==KAON && pid2==PION) ) { // Both the first and second are pion.
            // Combine two pions to caculate Ks mass
            float px = p1->px()+p2->px(); // x component of momentum
            float py = p1->py()+p2->py(); // y component of momentum
            float pz = p1->pz()+p2->pz(); // z component of momentum
            float e = p1->e()+p2->e();    // Energy
            float mass = sqrt ( e*e-px*px-py*py-pz*pz ); // Calculate the combined mass
            h_Dmass->Fill ( mass );  // Store mass in a histogram
            h_D_LorFit->Fill ( mass );
          }
        }
      }
    }
    // Clear at the end of particle loop.
    event->Clear();
  }
  // Let us show how many events are processed.
  printf( "***** %d events processed. Exit.\n", nevent_process );

  double bin1 = h_D_LorFit->FindFirstBinAbove(h_D_LorFit->GetMaximum()*3/4);
  double bin2 = h_D_LorFit->FindLastBinAbove(h_D_LorFit->GetMaximum()*3/4);
  double gamma = h_D_LorFit->GetBinCenter(bin2) - h_D_LorFit->GetBinCenter(bin1);
  double x0 = h_D_LorFit->GetBinCenter(h_D_LorFit->GetMaximumBin());
  double max = h_D_LorFit->GetMaximum();

  TF1 *L = new TF1("L",fitFunction,1.7,2.0,6);
  L->SetParameters(1,1,1,max,gamma,x0);
  h_D_LorFit->Fit("L");
//  h_Dmass->Draw("E");
  h_D_LorFit->Draw("E");
  gPad->Update();
  TPaveStats *st = (TPaveStats*)h_D_LorFit->FindObject("stats");
  gStyle->SetOptFit(0100);

  // Let us store  all the histogram data in the file named histfile.
  hf->Write();
}

// The function showhist is used when the histograms are displayed.
void showhist() {
  TFile *f = new TFile ( "histfile" );
  TBrowser *t = new TBrowser;
}
