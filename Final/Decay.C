// Author : Debajyoti Sengupa, Dept. of Physics, University of Calcutta


// How to Execute: 
// 1) Ensure you have zdecay.h file in the same directory as this.
// 2) Open a root session and type : .L Decay.C  
// 3) You should get a message stating that zdecay is being processed and roofit has been loaded
// 4) type : Decay t; (You may use any other var other than t)
// 5) t.Loop()  (This will execute your macro)

// The following code scans a tree in a root file and generates the z -> e-e+ peak from the data. 
// The code fits the signal and background regions using inbuilt roofit functions and displays the fit
// on a canvas.

// There a particle class defined at the very begining with the necessary accessor and mutator functions
// and the functionalities to convert the raw data into useful quantities. For example: Pseudo rapidity to theta.

// The compare_trans function is a boolean function that checks the transverse momentum (Pt) of two particles and returns true 
// if the first one is greater than the second one. Using this, the sort function sorts the vector of electrons according 
// to their Pt.

// RooFit PDFs are defined in the usual way. Upon execution of code, if you find HESSE does not converge, you may have to
// change the parameters a little. The fit will be displayed nonetheless. 
//Special thanks to Dr. Satyaki Bhattacharya and the Root forum for their help in writing this code.
//The code is also available at : https://github.com/DebajyotiS/SINP_ROOT


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


//* This is the main loop. This is what executes when you do .Loop() in the root session. Thus you put all your
//  macro code here

void Decay::Loop()
{

   
   cout << "Begin" << endl;
   if (fChain == 0) return;  
   
   TH1F *h1 = new TH1F("h1","e-e+ spectra;Mass;Events",200.,20.,500.);
   
   
   
   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"Total entries "<<nentries<<endl<<endl;
   Long64_t nbytes = 0, nb = 0; 
  

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      //if (ientry < 0) break;
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


          

         RooRealVar x("x","mass",20.,500.);
         RooDataHist dh("dh","e-e+",x,Import(*h1)) ;

         RooRealVar b("N_{backround}", "Number of background events", 0,4200);
         RooRealVar s("N_{signal}", "Number of signal events",0,700);

         x.setRange("Range1",65.,125.) ;
         x.setRange("Range2",20.,70.);
         x.setRange("Range3",130.,400.);
         
         
         // Defining the bernstein polynomial
         RooRealVar a1("a1","a1",-10,200);
         RooRealVar a2("a2","a2",-10,200);
         RooRealVar a3("a3","a3",-10,400);
         RooRealVar a4("a4","a4",-10,200);
         RooRealVar a5("a5","a5",-10,200);
         RooRealVar a6("a6","a6",-30,700);
         RooRealVar a7("a7","a7",-10,700);
         RooRealVar a8("a8","a8",-10,100);
         RooRealVar a9("a9","a9",-30,100);

         RooBernstein bg_bern("bg_bern","background",x, RooArgList(a1,a2,a3,a4,a5,a6,a7,a8,a9));
       

         // PDFs for convolution and CPDF
        
         //Defining the crystal ball PDF
         RooRealVar cmean("cmean","Central value of CB",0.,-30,30);
         RooRealVar csigma("#sigma","Width of CB",20,0,50.);
         RooRealVar calpha("#alpha","Alpha",7.,0,10);
         RooRealVar cn("cn","Order",6,-10.,20.);

         RooCBShape c_crystalball("c_crystalball", "convolution signal Region: CB", x,cmean,csigma,calpha,cn);


         //The Breit Wigner is very well known. So we just define the PDF explicitly. RooFit will see that it's parameters
         //are fixed.
         RooRealVar cmass("M_{Z}","Mass value for BW", 91.1876);
         RooRealVar cwidth("#sigma_{BW}","Spread of the BW", 2.49);

         RooBreitWigner c_BW("c_BW","convolution signal Region:BW",x,cmass, cwidth );


         RooFFTConvPdf CxB("CxB","CB (X) BW",x,c_crystalball,c_BW) ;
        

   
         RooAddPdf fullModel("fullModel", "CxB + bg_bern", RooArgList(CxB, bg_bern), RooArgList(s, b));

         RooFitResult* r = fullModel.fitTo(dh,Save()) ;
       
         
         RooPlot* frame1 = x.frame(Title("Imported Histogram and fullModel"));
         RooPlot* frame2 = x.frame(Title("PDF: fullModel = CB X Breit Wigner + Bernstein"));
         RooPlot* frame3 = x.frame(Title("Signal Only + Convolution fit"));
         RooPlot* frame4 = x.frame((Title("Background Model")));

         dh.plotOn(frame1,DataError(RooAbsData::None));
         fullModel.plotOn(frame1,LineColor(kRed));

         fullModel.plotOn(frame2, LineColor(kRed),Normalization(b.getVal()+s.getVal()));
         bg_bern.plotOn(frame2,LineColor(kBlack),LineStyle(2),Normalization(b.getVal()));

                 
          
         TH1 *h2 = bg_bern.createHistogram("background;mass;count", x, Binning(200,20.,500.));
      
         h1 -> Add(h2, -1*b.getVal());
         h1->SetMinimum(15);

         RooDataHist cdh("cdh","e-e+ signal",x,Import(*h1)) ;

         RooFitResult *rr = CxB.fitTo(cdh, Range("Range1"),Save());

         cdh.plotOn(frame3,DataError(RooAbsData::None));
         CxB.plotOn(frame3,LineColor(kRed));
         CxB.paramOn(frame3);

         dh.plotOn(frame4,DataError(RooAbsData::None));
         bg_bern.plotOn(frame4,LineColor(kGreen),LineStyle(2));
         
         TCanvas* c = new TCanvas("canvas","canvas",900,720) ;

         
         
         c -> Divide(2,2);
         c -> cd(1);
         frame1 -> Draw();
         c -> cd(2);
         frame2->Draw();
         c -> cd(3);
         frame4->Draw();
         c -> cd(4);
         frame3->Draw();
         

         // TFile f("result.root","recreate");
         // c -> Write();
         // r -> Write();
         // f.Close();
}









      // Old fit functions and saving routine

      //! The crystalbal+expo fit routine
      // TF1 *fit = new TF1("fit", "crystalball(0) + expo(5)", 0., 1000.);
      // fit ->SetParameters(130.0,90.0,9.0,1.0,300.0, 0.01);
      // for (Int_t i = 0; i < fit->GetNpar(); i++) std::cout << i << " : " << fit->GetParName(i) << std::endl;
      // h1->Fit("fit", "", "", 30., 200.);
      

      // TF1 *fline = new TF1("fline","expo(0)",30.,550.);
      // fline->SetParameters(300.,-0.01);
   
      // h1->Fit("fline","","",30.,65.);
      // h1->Fit("fline","","",105.,300.);
      

      //! This part is to create a file and save the plots.

      // TFile f("Histo.root","recreate");
      // h1 -> Write();
      // fit -> Write();
      // f.Close();
      // c1 -> cd();
      // h1 ->Draw();
      





// Roofit with expo and gaussian

//  RooRealVar x("x","mass",0.,1000.);
//    // x.setRange("Range1",35.75,50.13) ;
//    x.setRange("Range2",110.101,200.) ;
//    RooRealVar lambda("lambda", "slope", -0.001, -1., 1.);
//    RooExponential expo("expo", "expo", x, lambda);
//    RooExponential expo1("expo1", "expo", x, lambda);
//    RooDataHist dh("dh","e-e+ peak histo",x,Import(*h1)) ;
//    // RooFitResult* r = 
//    RooExtendPdf background("background","Side Bands",expo,x,"FULL") ;
//    RooFitResult* r = expo.fitTo(dh,Range("Range2"), Save(kTRUE)) ;
//    r -> Print();
   
//    RooPlot* frame = x.frame() ;
//    RooPlot* frame1 = x.frame() ;

      
//    dh.plotOn(frame,DrawOption("B"),DataError(RooAbsData::None),FillColor(kGray));
//    expo.plotOn(frame);
//    // background.plotOn(frame,LineStyle(kDashed),LineColor(kRed));
   

//    RooDataSet *h2 = expo.generate(x,2500);
   
   
//    h2-> plotOn(frame1, LineColor(kBlue));
   
//    TH1* h2_hist = h2->createHistogram("h2_hist", x);
//    RooDataHist dh2("dh2","e-e+ peak histo",x,Import(*h2_hist)) ;
//    RooFitResult* r2 = expo1.fitTo(dh2, Save(kTRUE)) ;
//    r2 -> Print();
   
//    RooPlot *frame3 = x.frame();
//    dh2.plotOn(frame3);
//    expo1.plotOn(frame3);
   
//    TCanvas* c = new TCanvas("canvas","canvas",800,800) ;
//    c -> Divide(2,2);
//    c -> cd(1);
//    frame -> Draw();
//    c -> cd(2);
//    frame1 -> Draw();
//    c -> cd(3);
//    frame3 -> Draw();


// CB (Gaus) + Bernstein Polynomial fit :
// RooRealVar x("x","mass",30.,500.);
//x.setRange("Range1",25.,75.) ;
//x.setRange("Range2",107.,500.) ;
//RooRealVar a1("a1","a1",-50,1000);
//RooRealVar a2("a2","a2",-50,400);
//RooRealVar a3("a3","a3",-50,400);
//RooRealVar a4("a4","a4",-500,1000);
//RooRealVar a5("a5","a5",-50,1000);
//RooRealVar a6("a6","a6",-50,1000);
//RooRealVar a7("a7","a7",-50,1000);
//RooRealVar a8("a8","a8",-50,100);
//RooRealVar a9("a9","a9",-50,100);
//RooRealVar aa("aa","aa",-50,100);

//RooDataHist dh("dh","e-e+ peak histo",x,Import(*h1)) ;

//RooRealVar lambda("lambda", "slope", -0.01, -1., 1.);
//RooExponential expo("expo", "expo", x, lambda);
  

//RooRealVar mass("mass","Central value of Gaussian",90.,80.,100.);
//RooRealVar sigma("sigma","Width of Gaussian",20,0,100);
//RooRealVar alpha("alpha","Alpha",40.,0.,100.);
//RooRealVar n("n","Order",6,0.,10.);

//RooCBShape crystalball("crystalball", "Signal Region", x,mass,sigma,alpha,n);

//RooGaussian signal("gaus", "The signal distribution", x, mass, sigma); 

//RooRealVar b("b", "Number of background events", 0, 5500);
//RooRealVar s("s", "Number of signal events", 0, 500);

//RooBernstein bg_bern("bg_bern","background",x, RooArgList(a1,a2,a3,a4,a5,a6,a7,a8,a9));

//// RooAddPdf fullModel("fullModel", "signal + bg_bern", RooArgList(signal, bg_bern), RooArgList(s, b));
//RooAddPdf fullModel("fullModel", "crystalball + bg_bern", RooArgList(crystalball, bg_bern), RooArgList(s, b));

//RooFitResult* r = fullModel.fitTo(dh,Save()) ;
//r -> Print();
  
//RooPlot* frame1 = x.frame(Title("Imported Histogram and fullModel fit"));
//RooPlot* frame2 = x.frame(Title("PDF: fullModel = CB + Bernstein"));
//RooPlot* frame3 = x.frame(Title("Background: Bernstein Polynomial (N = 8)"));

//dh.plotOn(frame1,DataError(RooAbsData::None));
//// bg_bern.plotOn(frame1);
//fullModel.plotOn(frame1,LineColor(kRed));
//fullModel.plotOn(frame2,LineColor(kRed));
//bg_bern.plotOn(frame3,LineColor(kBlack));
//// expo.plotOn(frame);

//fullModel.paramOn(frame2,Layout(0.55,0.95,0.8)) ;
  
//// frame -> Draw();
//TCanvas* c = new TCanvas("canvas","canvas",1366,768) ;
//c -> Divide(2,2);
//c -> cd(1);
//frame1 -> Draw();
//c -> cd(2);
//frame2 ->Draw();
//c -> cd(3);
//frame3 ->Draw();
//c -> cd(4);
//h1 -> Draw();

      
      
     
   
