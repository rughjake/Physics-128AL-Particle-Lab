// File : BEvent.cc
// Description : Implementation of BEvent class // 
// Author : Ryosuke Itoh, IPNS, KEK // Date : 28 - Jan - 2004
//-

 #include "BParticle.h"
 #include "BEvent.h"

 ClassImp(BEvent)

 BEvent::BEvent()
 {
   m_particles = new TClonesArray ( "BParticle", 500 );
   m_evno = 0;
 }

 BEvent::~BEvent()
 {
   delete m_particles;
 }

 void BEvent::EventNo ( int evt )
 {
   m_evno = evt;
 }

 int BEvent::EventNo()
 {
   return m_evno;
 }

 void BEvent::AddTrack ( float px, float py, float pz, float e,
 			float charge, SIMPLEPID pid )
 {
   TClonesArray& particles = *m_particles;
   new (particles[m_nprt++]) BParticle ( px, py, pz, e, charge, pid ); 
 }

 int BEvent::NParticles()
 {
   return m_nprt;
 }

 TClonesArray* BEvent::GetParticleList () {
   return m_particles;
 }

 void BEvent::Clear ( )
 {
   m_particles->Clear();
   m_nprt = 0;
 }

