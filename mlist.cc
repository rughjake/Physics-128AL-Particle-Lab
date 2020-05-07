// The numbers of particles in each event are stored in a histogram.
// All particle masses in the analyzed events are calculated and stored in a histogram.
//+
// File : mlist.cc
// Description : Analyze a data file
//
// Author : Ryosuke Itoh, IPNS, KEK
// Date : 1 - Mar - 2004
//-
// Modified by T.Nozaki  (2005/10/04)
// Modified by S.Nishida (2007/02/03)
// Modified by T.Nozaki (2009/07/01) Comments in English

// Modified by Jake Rugh (5/5/2020) Debugged, added histogram axis label
#include <cstdio>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TClonesArray.h"
#include "TBrowser.h"

#include "BParticle.h"
#include "BEvent.h"
//enum SIMPLEPID {PHOTON, ELECTRON, PION, MUON, KAON, PROTON};

using namespace std;

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
  h_nprt->GetXaxis()->SetTitle("particle count per event");
  // The distribution of all particle masses.
  // The longitudinal axis is divided into 1000 bins between 0.0 and 1.0 GeV 
  TH1F *h_mass1 = new TH1F ( "h_mass1", "mass of all particles",
			     1000, 0.0, 1.0 );
  h_mass1->GetXaxis()->SetTitle("mass (GeV)");
  // The longitudinal axis is divided into 1000 bins between 0.0 and 0.001 GeV 
  TH1F *h_mass2 = new TH1F ( "h_mass2", "mass of all particles",
			     1000, 0.0, 0.001);
  h_mass2->GetXaxis()->SetTitle("mass (GeV)");

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
 
    // Let us extract the number of particles in each event.
    int nparticles = event->NParticles();
    // It is stored in the histogram named h_nprt.
    h_nprt->Fill ( nparticles );

    // Particles are treated as an array.
    TClonesArray& plist = *(event->GetParticleList());

    // loop of all particles.
    for ( int j=0;j<nparticles;j++ ) {
      // The j-th particle is named  p1.
      BParticle* p1 = (BParticle*)plist[j];
      // Obtain the momentum and energy of a particle.
      float px = p1->px();  // x component of momentum.
      float py = p1->py();  // y component of momentum.
      float pz = p1->pz();  // z component of momentum.
      float  e = p1->e();   // Energy.
      float mass = sqrt ( e*e-px*px-py*py-pz*pz ); // Calculate the mass.
      h_mass1->Fill ( mass );  // Store mass in the histogram named h_mass1.
      h_mass2->Fill ( mass );  // Store mass in the histogram named h_mass2.
    }
    // Clear at the end of particle loop.
    event->Clear();
  }
  // Let us show how many events are processed.
  printf( "***** %d events processed. Exit.\n", nevent_process );
  // Let us store  all the histogram data in the file named histfile.

  hf->Write();
 }

// The function showhist is used when the histograms are displayed.
void showhist() {
  TFile *f = new TFile ( "histfile" );
  TBrowser *t = new TBrowser;
}
