#ifndef BPARTICLE_H
#define BPARTICLE_H
//+
// File : BParticle.h
// Description : class to contain particle info.
//
// Author : Ryosuke Itoh, IPNS, KEK
// Date : 28 - Jan - 2004
//-

#include "TObject.h"

enum SIMPLEPID {PHOTON, ELECTRON, PION, MUON, KAON, PROTON};

class BParticle : public TObject {

 private:
  float m_px;
  float m_py;
  float m_pz;
  float m_e;
  float m_charge;
  SIMPLEPID m_pid;

 public:
  BParticle();
  BParticle ( float px, float py, float pz, float e, float c, SIMPLEPID pid );
  ~BParticle();

  float px();
  float py();
  float pz();
  float e();
  float charge();
  SIMPLEPID pid();

  ClassDef ( BParticle, 1 ) // Simple particle class
};
#endif

