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
   cout << "Begin" << endl;
   if (fChain == 0) return;  
   
   TH1F *h1 = new TH1F("h1","e-e+ spectra;Mass;Events",200,30.,500.);
   
   
   
   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"Total entries "<<nentries<<endl<<endl;
   Long64_t nbytes = 0, nb = 0; 
  

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      
      //* Creating a vector of objects of class particles.
      vector<particles> electrons;
  
      
      // Acquisition loop. If the event has at least one electron (or particle) then it makes an electron object
      // which is inserted into the vector - electrons

      if(nEle>1) {
         
         for (int i=0 ;i<nEle; i++){
            particles electron(elePhi -> at(i),eleEta->at(i),elePt->at(i),eleCharge->at(i));
            electrons.push_back(electron);
         						   }
         std::sort(electrons.begin(),electrons.end(), compare_trans);
         
         // Contructing the four vector for the first two particles with highest Pt - Summing them in V3 for display.
         TLorentzVector v1;
         v1.SetXYZT(electrons[0].getPx(),electrons[0].getPy(),electrons[0].getPz(),electrons[0].getE());
         TLorentzVector v2;
         v2.SetXYZT(electrons[1].getPx(),electrons[1].getPy(),electrons[1].getPz(),electrons[1].getE());
         TLorentzVector v3 = v1+v2;
         h1 -> Fill(v3.M());
      
                 }
      
                                                    }
      
 }









