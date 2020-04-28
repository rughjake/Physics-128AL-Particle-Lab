// File : 128.cc
// Author : Jake Rugh
//
// Analyze all combinations of decay particles to search for presence of parent particles
//
// Edited from :
//
// File : mlist.cc
// Author : Ryosuke Itoh, IPNS, KEK
// Modified by T.Nozaki  (2005/10/04)
// Modified by S.Nishida (2007/02/03)
// Modified by T.Nozaki (2009/07/01) Comments in English

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

  TFile *f = new TFile( file );
  TTree *t = (TTree*)f->Get("T");
  BEvent* event = new BEvent();
  TBranch *branch = t->GetBranch("BEvent");
  branch->SetAddress(&event);

  // Initialization codes finish here.


  // New histogram folder. The first quoted string in the parantheses is the name.
  TFile *hf = new TFile("histfile", "RECREATE" );

  /* New histogram file. First argument is file name, second argument is descriptor,
     third argument is number of bins, fourth and fifth arguments are bin value range. */
  TH1F *Hist1 = new TH1F ( "Hist1", "mass of all hypothetical parent particles", 1000, 0, 1 );


  // Number of events in data
  int nevent = (int)t->GetEntries();

  // Only analyzes certain number of events if specific by user.
  // Else, analyzes all events
  int nevent_process = nevent;
  if( maxevt>0 && maxevt<nevent )
    nevent_process = maxevt;

  // loop of events.
  for ( int i=0;i<nevent_process;i++ ) {
    //  write an event number every 1000 events.
    if ( i%1000 == 0 ) printf ( "***** Event %d\n", i );

    // An event is read.
    int nb = t->GetEntry ( i );

    // number of particles in each event.
    int nparticles = event->NParticles();

    // Particles entered into array.
    TClonesArray& plist = *(event->GetParticleList());

    // Analyzing all combinations of decay particles for each event 
    for ( int m=0;m<nparticles;m++ ) {
      BParticle* p1 = (BParticle*)plist[m];
        float px1 = p1->px(); // x momentum
        float py1 = p1->py(); // y momentum
        float pz1 = p1->pz(); // z momentum
        float  e1 = p1->e();  // energy
      for ( int n=0;n<nparticles;n++ ) {
        if ( m != n) // I do not want to analyze particles against themselves
          p1 = (BParticle*)plist[n];
          float px2 = p1->px();  // x momentum
          float py2 = p1->py();  // y momentum
          float pz2 = p1->pz();  // z momentum
          float  e2 = p1->e();   // energy

          // calculate momentum and energy of hypothetical parent particle,
          // using conservation of momentum and energy
          float px3 = px1 + px2;
          float py3 = py1 + py2;
          float pz3 = pz1 + pz2;
          float  e3 = e1 + e2;
          float mass = sqrt ( e3*e3-px3*px3-py3*py3-pz3*pz3 ); // mass of parent
          Hist1->Fill( mass ); // enter parent mass into histograms
        }
    }
    // clear event
    event->Clear();
  }
  // print number of events processed
  printf( "***** %d events processed. Exit.\n", nevent_process );
  // store histogram data in histfile
  hf->Write();
  // display histograms
  Hist1->Draw();
 }

// The function showhist is used when the histograms are displayed.
void showhist() {
  TFile *f = new TFile ( "histfile" );
  TBrowser *t = new TBrowser;
}
