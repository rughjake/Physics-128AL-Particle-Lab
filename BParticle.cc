//+
// File : BParticle.cc
// Description : Implementation of BTrack class //
// Description : Ryosuke Itoh, IPNS, KEK // Date : 28 - Jan - 2004
//-

#include "BParticle.h"

ClassImp (BParticle)

BParticle::BParticle()
{
 }

 BParticle::BParticle ( float px, float py, float pz, float e,
 		       float charge, SIMPLEPID pid ) {
   m_px = px;
   m_py = py;
   m_pz = pz;
   m_e = e;
   m_charge = charge;
   m_pid = pid;
 }

 BParticle::~BParticle()
 {
 }

 float BParticle::px ()
 {
   return m_px;
 }
 float BParticle::py ()
 {
   return m_py;
 }
 float BParticle::pz ()
 {
   return m_pz;
 }
 float BParticle::e ()
 {
   return m_e;
 }
 float BParticle::charge ()
 {
   return m_charge;
 }
 SIMPLEPID BParticle::pid ()
 {
   return m_pid;
 }
