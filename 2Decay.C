#define Decay_cxx
#include "Decay.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TFormula.h>
#include <TF1.h>
#include <TList.h>

// Writing a Particle Class

class particles{
public:

   //*Constructor
   particles(double elePhi, double eleEta, double elePt, double eleCharge){
      phi = elePhi;
      eta = eleEta;
      Pt = elePt;
      charge = eleCharge;
   }

   //* destructor
   ~particles(){
    }

    //* Accessors
   double getphi(){
      return phi;
   }
   double geteta(){
      return eta;
   }
   double getPt(){
      return Pt;
   }
   double getE(){
      return Pt/(sin(2.0*atan(exp(-eta))));
   }
   double getPx(){
      return Pt*cos(phi);
   }
   double getPy(){
      return Pt*sin(phi);
   }
   double getPz(){
      return Pt*cos(phi)/sin(phi);
   }
   double getcharge(){
      return charge;
   }


   //* Mutators
   void setphi(double Phi){
      phi = Phi;
   }
   void seteta(double Eta){
      eta = Eta;
   }
   void setPt(double Pt){
      Pt = Pt;
   }
   void setcharge(double Charge){
      charge = Charge;
   }

private:
   double E,phi,theta,eta,Pt,Px,Py,Pz,charge;
};


//* This is the comparing function. It takes two values and returns False/True based on the expression used.
bool compare_trans(particles &particle1, particles &particle2){
    return particle1.getPt() > particle2.getPt();
   }




void Decay::Loop()
{

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}
