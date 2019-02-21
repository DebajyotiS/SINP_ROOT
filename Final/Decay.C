#define Decay_cxx
#include "Decay.h"
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TFormula.h>
#include <TF1.h>
#include <TList.h>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include <RooBernstein.h>
#include <RooBreitWigner.h>
#include <RooCBShape.h>
using namespace RooFit ;


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



    //Defining the observables and parameters and initialising the values                                               
    RooRealVar x("x","mass",30.,500.);
	RooRealVar mass("mass","Central value of CB",90.18,90.17,90.19);
    RooRealVar sigma("sigma","Width of CB",20,0,100);
    RooRealVar alpha("alpha","Alpha",40.,-100,100.);
    RooRealVar n("n","Order",6,0.,10.);

    //Defining the crystalball pdf
    RooCBShape crystalball("crystalball", "Signal Region", x,mass,sigma,alpha,n);    


    RooRealVar a1("a1","a1",-50,1000);
    RooRealVar a2("a2","a2",-50,400);
    RooRealVar a3("a3","a3",-50,400);
    RooRealVar a4("a4","a4",-500,1000);
    RooRealVar a5("a5","a5",-50,1000);
    RooRealVar a6("a6","a6",-50,1000);
    RooRealVar a7("a7","a7",-50,1000);
    RooRealVar a8("a8","a8",-50,100);
    RooRealVar a9("a9","a9",-50,100);

    //Defining the bernstein polynomial pdf.
	RooBernstein bg_bern("bg_bern","background",x, RooArgList(a1,a2,a3,a4,a5,a6,a7,a8,a9));   

	//Defining the number of background and signal events. These will the weights with which your pdfs get added.
	RooRealVar b("b", "Number of background events", 0, 5500);
    RooRealVar s("s", "Number of signal events", 0, 600); 

    //Finally this is the "fullmodel pdf"
    RooAddPdf fullModel("fullModel", "crystalball + bg_bern", RooArgList(crystalball, bg_bern), RooArgList(s, b));

    //Import the TH1F histogram as a RooDataHist
    RooDataHist dh("dh","e-e+",x,Import(*h1)) ;

    //Fit the histogram with the fullmodel pdf 
    RooFitResult* r = fullModel.fitTo(dh,Save()) ;







    RooPlot* frame1 = x.frame(Title("Imported Histogram and fullModel fit"));
    RooPlot* frame2 = x.frame(Title("PDF: fullModel = CB + Bernstein"));
    RooPlot* frame3 = x.frame(Title("Background: Bernstein Polynomial (N = 8)"));
    RooPlot* frame4 = x.frame(Title("Background: Histogram"));

    dh.plotOn(frame1,DataError(RooAbsData::None));
    
    fullModel.plotOn(frame1,LineColor(kRed));
    fullModel.plotOn(frame2,LineColor(kRed));
    bg_bern.plotOn(frame2,LineColor(kGreen));
	fullModel.paramOn(frame1,Layout(0.55,0.95,0.8));

	TCanvas* c = new TCanvas("canvas","canvas",1024,1024) ;
	c -> Divide(2,2);
    c -> cd(1);
    frame1 -> Draw();
    c -> cd(2);
    frame2->Draw();
    c -> cd(3);
    h1 -> Draw();
    c -> cd(4);
    frame3->Draw();



      
 }









