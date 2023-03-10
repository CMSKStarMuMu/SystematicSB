//
//
// 
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sys/stat.h>
#include "Riostream.h"
#include <map>
#include <string>
//#include <boost/algorithm/string.hpp>
//#include <boost/algorithm/string/trim.hpp>
#include <vector>
#include <math.h>
//#include <TCint.h>
//#include <TGenericClassInfo.h> 
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include <TROOT.h>
#include <TEnv.h>
#include <TSystem.h>
#include <TTree.h>
#include "TBranch.h"
#include <TApplication.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h> 
#include <TF1.h>  
#include <TF2.h> 
#include <TLorentzVector.h>
#include <TVector3.h>
#include "TChain.h"
#include <time.h> 
#include <TSystemDirectory.h>
#include <TEfficiency.h>
#include <TLegend.h>
#include <TMinuit.h>
#include "Math/WrappedMultiTF1.h"
#include "TRandom.h" 
#include "TRandom3.h" 
#include  <TStopwatch.h>
#include "TH1F.h"
#include "TH2F.h"			// unused?
#include "TStyle.h"
#include "TCanvas.h"
#include <TGraphAsymmErrors.h>
#include <TFrame.h>
#include <TFitResult.h>
#include <TFitter.h>
#include "Fit/Fitter.h"
#include <TMatrixDSym.h>
#include <TMatrixD.h>
#include <TBinomialEfficiencyFitter.h>
#include <TKDTreeBinning.h>
#include <TH2Poly.h>
//#if !defined(__CINT__) || defined(__MAKECINT__)
#include <RooRandom.h>
#include <RooFit.h>
#include <RooMinuit.h>
#include <RooRealVar.h>
#include <RooConstVar.h>
#include <RooPlot.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooRealVar.h>
#include <RooAbsReal.h>
#include <RooArgList.h>
#include <RooRealProxy.h>
#include <RooCategoryProxy.h>
#include <RooAbsCategory.h>
#include <RooProdPdf.h>
#include <RooAddPdf.h>
#include <RooProduct.h>
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooPolynomial.h>
#include <RooChebychev.h>
#include <RooWorkspace.h>
#include <RooExponential.h>
#include <RooErrorVar.h>
#include <RooFitResult.h>
#include <RooRangeBinning.h>
#include <RooBinning.h>
#include <RooNumGenConfig.h>
#include <RooBernstein.h>
#include <RooPolynomial.h>
#include <RooExtendPdf.h>
#include <RooSimultaneous.h>
#include <RooMultiVarGaussian.h>
#include "RooBernsteinSideband.h"
#include "RooMCStudy.h"
#include "TRatioPlot.h"

#include <sys/time.h>
#include <sys/times.h>
#include <iostream>
// System stuff
#include <fstream> 
#include <sys/time.h>
#include <sys/times.h>

timeval startTime, stopTime, totalTime;
clock_t startCPU, stopCPU; 
tms startProc, stopProc; 

using namespace std; 
using namespace ROOT;
using namespace RooFit;
void ReadSimultaneousModel();
char OutSaveFileName[400] =  "/gwpool/users/dini/Fittone/UML-fit-MT-bin6-GenSB/simFitResults4d/simFitResult_data_fullAngularMass_Swave_2016_MCStat_b6_fullStat_noMeanCon.root";
char *ReadOutSaveFileName;

//std::string PNGNameReadSB3D = "systematicSideband" ;

int    xCosLHBin =  25;
int    xCosKHBin =  25;
int    xPhiHBin  =  25;
double XMinCosThetaL	     = -1.;
double XMaxCosThetaL	     =  1.;
double XMinCosThetaK	     = -1.;
double XMaxCosThetaK	     =  1.;
double XMinPhi  	     =-TMath::Pi();
double XMaxPhi  	     = TMath::Pi();
double XMinSign = 5.0;
double XMaxSign = 5.6;
double XStepSign = 0.025;
float  xMassHBin = (XMaxSign -XMinSign)/XStepSign;

int Maxsystematic = 100;

int Q2Bin  =0 ;
int RunEra =0 ;

int GenEntries=1000;
//=================================================
int main (int argc, char** argv) {
  TStopwatch TimeWatch;
  TimeWatch.Start();
  
  if (argc>1 ){
     ReadOutSaveFileName  = (char *) malloc(strlen(argv[1])+1);
     strcpy(ReadOutSaveFileName,argv[1]);
     sprintf(OutSaveFileName,"%s",ReadOutSaveFileName);
     std::cout<<"========================================================================="<<endl;
     std::cout<<Form("Setting the File: %s",argv[1])<<std::endl;
     std::cout<<"========================================================================="<<endl;
  }else{
     std::cout<<"========================================================================="<<endl;
     std::cout<<Form("Setting DEFAULT File: %s",OutSaveFileName)<<std::endl;
     std::cout<<"========================================================================="<<endl;
  }
//  PNGNameReadSB3D=Form("%s-%s.png",PNGNameReadSB3D.c_str(),OutSaveFileName);
  

  ReadSimultaneousModel(); 
//  app.Run() ;
  cout<<"esco..." <<endl;
  TimeWatch.Stop();
  TimeWatch.Print();
  
  return 0 ;
}
void  ReadSimultaneousModel(){
    RooRealVar* ctL = new RooRealVar("ctL", "ctL",  XMinCosThetaK,XMaxCosThetaK);
    RooRealVar* ctK = new RooRealVar("ctK", "ctK",  XMinCosThetaL,XMaxCosThetaL);
    RooRealVar* phi = new RooRealVar("phi", "phi",  XMinPhi,XMaxPhi);
    RooArgList observables (*ctK, *ctL, *phi);
//    RooRealVar* tagged_mass = new RooRealVar("tagged_mass" , "#mu^{+}#mu^{-}K#pi mass", XMinSign,XMaxSign, "GeV/c^2");
    
    TFile* OutSaveFile = new TFile( OutSaveFileName, "READ" );
    if ( !OutSaveFile || !OutSaveFile ->IsOpen() ) {
      cout<<Form("File not found: %s\n",OutSaveFileName)<<endl;
      exit(1);
    }
    RooWorkspace* w = (RooWorkspace*)OutSaveFile->Get("wsp_out");
    if ( !w || w->IsZombie() ) {
     cout<<Form("Workspace not found in file:%s\n",OutSaveFileName)<<endl;
     exit(1);
    } else {
     cout<<Form("Workspace Found!!! In file : %s\n",OutSaveFileName)<<endl;
    }
    RooSimultaneous* simPdf = (RooSimultaneous*)w->pdf("simPdf");
    RooArgSet *params      = (RooArgSet *)simPdf->getParameters(observables);
    auto iter  = params->createIterator();
    RooRealVar* ivar =  (RooRealVar*)iter->Next();
    while (ivar) {
      cout <<Form("%s=%f +/-%f",ivar->GetName(),ivar ->getVal(),ivar ->getError())<<endl;
      ivar =  (RooRealVar*)iter->Next();
    } 
  } 
