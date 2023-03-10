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
void ReadSBModel();
char OutSaveFileName[400] =  "savesb-2018-Q2Bin-3-Bins-25-25-25-BernDeg-5-4-3-WSBL-5-2dot2-WSBR-2dot4-5dot6-SigmaProb.root";
char *ReadOutSaveFileName;

char PNGNameReadSB3D[400] = "renameSideband";
RooDataSet *RenameeDataNew=0;
RooDataSet *RenameeMassNew=0;
RooBernsteinSideband *bkg_ang_pdf=0;
RooAbsPdf *bkg_mass_sb=0;
//RooExponential *bkg_mass_sb=0;
RooRealVar *max_sbl=0;
RooRealVar *min_sbr=0;
TMatrixD *covMatrix =0;

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

int Maxrename = 100;

int Q2Bin  =0 ;
int RunEra =0 ;

int RenameEntries=1000;
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
  sprintf(PNGNameReadSB3D,"%s-%s.png",PNGNameReadSB3D,OutSaveFileName);
  

  ReadSBModel(); 
//  app.Run() ;
  cout<<"esco..." <<endl;
  TimeWatch.Stop();
  TimeWatch.Print();
  
  return 0 ;
}
void  ReadSBModel(){
    RooRealVar* ctL = new RooRealVar("ctL", "ctL",  XMinCosThetaK,XMaxCosThetaK);
    RooRealVar* ctK = new RooRealVar("ctK", "ctK",  XMinCosThetaL,XMaxCosThetaL);
    RooRealVar* phi = new RooRealVar("phi", "phi",  XMinPhi,XMaxPhi);
//    RooRealVar* tagged_mass = new RooRealVar("tagged_mass" , "#mu^{+}#mu^{-}K#pi mass", XMinSign,XMaxSign, "GeV/c^2");
    
    TFile* OutSaveFile = new TFile( OutSaveFileName, "READ" );
    if ( !OutSaveFile || !OutSaveFile ->IsOpen() ) {
      cout<<Form("File not found: %s\n",OutSaveFileName)<<endl;
      exit(1);
    }
    RooWorkspace* w = (RooWorkspace*)OutSaveFile->Get("wsb");
    if ( !w || w->IsZombie() ) {
     cout<<Form("Workspace not found in file:%s\n",OutSaveFileName)<<endl;
     exit(1);
    } else {
     cout<<Form("Workspace Found!!! In file : %s\n",OutSaveFileName)<<endl;
    }
//    w->Print();
    
//
//  Cov Matrix
//   
//     TMatrixD *covMatrix = (TMatrixD*)w->obj(Form("covMatrix_bin%d_%d",Q2Bin,RunEra));
//     if(covMatrix==0){
//         cout<<Form("covMatrix_bin%d_%d not found in file:%s\n",OutSaveFileName)<<endl;
//     }else{  
//      cout<<Form("covMatrix found in file:%s  DUMP:\n",OutSaveFileName)<<endl;
//      covMatrix->Print("f=  %10.3e  ");
//     } 
    bool bfound=false;
    for(int i=0;i<8;i++){
     if (bfound) break;
     for(int j=6;j<9;j++){
      if (bfound) break;
      covMatrix = (TMatrixD*)w->obj(Form("covMatrix_bin%d_201%d",i,j));
      if ( !covMatrix) {
       if(i==7&&j==8){
//         cout<<Form("covMatrix_bin%d_201%d not found in file:%s\n",i,j,OutSaveFileName)<<endl;
//         exit(1);
       }	 
      } else {
       cout<<Form("covMatrix_bin%d_201%d Found!! in file:%s DUMP:\n",i,j,OutSaveFileName)<<endl;
       bfound=true;
      }
     }  
    }
    if(bfound) {
     covMatrix->Print("f=  %10.3e  ");
    }else{
     cout<<Form("NO covMatrix found in file:%s !!!\n",OutSaveFileName)<<endl;
    } 
//  
    bfound=false;
    for(int i=0;i<8;i++){
     if (bfound) break;
     for(int j=6;j<9;j++){
      if (bfound) break;
      bkg_ang_pdf= (RooBernsteinSideband *)w->pdf(Form("BernSideBand_bin%d_201%d",i,j));
      if ( !bkg_ang_pdf ) {
//       cout<<Form("BernSideBand_201%d not found\n",i)<<endl;
       if(i==7&&j==8){
         cout<<Form("Not found any BernSideBand function in  in file:%s\n",OutSaveFileName)<<endl;
//         cout<<Form("BernSideBand_bin%d_201%d not found in file:%s\n",i,j,OutSaveFileName)<<endl;
         exit(1);
       }	 
      } else {
       cout<<Form("BernSideBand_bin%d_201%d Found!! in file:%s\n",i,j,OutSaveFileName)<<endl;
       Q2Bin=i;
       RunEra=2010+j;
       bfound=true;
      }
     }  
    }
    bfound=false;
    for(int i=0;i<=7;i++){
     if (bfound) break;
     for(int j=6;j<=8;j++){
      if (bfound) break;
      bkg_mass_sb = (RooAbsPdf *)w->pdf(Form("bkg_mass_sb_bin%d_201%d",i,j));
//      bkg_mass_sb = (RooExponential *)w->pdf(Form("bkg_mass_sb_bin%d_201%d",i,j));
      if ( ! bkg_mass_sb) {
       if(i==7&&j==8){
         cout<<Form("bkg_mass_sb_bin*_201* not found in file:%s\n",OutSaveFileName)<<endl;
         exit(1);
       }
      } else {
       cout<<Form("bkg_mass_sb_bin%d_201%d Found!! in file:%s\n",i,j,OutSaveFileName)<<endl;
       bfound=true;
      } 
     }
    }
    
    
    
    bfound=false;
    for(int i=0;i<=7;i++){
     if (bfound) break;
     for(int j=6;j<=8;j++){
      if (bfound) break;
      max_sbl = (RooRealVar *)w->var(Form("max_sbl_bin%d_201%d",i,j));
      if ( ! max_sbl) {
       if(i==7&&j==8){
         cout<<Form("max_sbl_bin*_201* not found in file:%s\n",OutSaveFileName)<<endl;
         exit(1);
       }
      } else {
       cout<<Form("max_sbl_bin%d_201%d=%f\n",i,j,max_sbl->getVal())<<endl;
       bfound=true;
      } 
     }
    }
    bfound=false;
    for(int i=0;i<=7;i++){
     if (bfound) break;
     for(int j=6;j<=8;j++){
      if (bfound) break;
      min_sbr = (RooRealVar *)w->var(Form("min_sbr_bin%d_201%d",i,j));
      if ( ! min_sbr ) {
       if(i==7&&j==8){
         cout<<Form("min_sbr_bin*_201* not found in file:%s\n",OutSaveFileName)<<endl;
         exit(1);
       }
      } else {
       cout<<Form("min_sbr_bin%d_201%d=%f\n",i,j,min_sbr->getVal())<<endl;
       bfound=true;
      } 
     }
    }

    int  maxDegree1 =   bkg_ang_pdf->maxDegree1Get();
    int  maxDegree2 =   bkg_ang_pdf->maxDegree2Get();
    int  maxDegree3 =   bkg_ang_pdf->maxDegree3Get();
    int  numParameters = (maxDegree1+1)*(maxDegree2+1)*(maxDegree3+1);

    cout<<Form("maxDegree1=%d\n",maxDegree1); 
    cout<<Form("maxDegree2=%d\n",maxDegree2); 
    cout<<Form("maxDegree3=%d\n",maxDegree3); 
    cout<<Form("numParameters=%d\n",numParameters); 

    RooArgList observables (*ctK, *ctL, *phi);
//     std::vector<TH1D*> HPar;
//     std::vector<TH3D*> HBern;
//     std::vector<TH1D*> HBernX;
//     std::vector<TH1D*> HBernY;
//     std::vector<TH1D*> HBernZ;
//     std::vector<RooBernsteinSideband*> BernSideBand; 
     RooArgList* bkg_ang_params = (RooArgList*)bkg_ang_pdf->getParameters(observables);
//    auto iter = bkg_ang_params->createIterator();
//    RooRealVar* ivar =  (RooRealVar*)iter->Next();
//     RooArgList* bkg_ang_params_rename = new RooArgList("bkg_ang_params_rename");
//     RooArgList* bkg_ang_params_redu_tmp = new RooArgList("bkg_ang_params_redu_tmp");
//     RooArgList* bkg_ang_params_rename_redu_tmp = new RooArgList("bkg_ang_params_rename_redu_tmp");
//     RooRealVar * rndPar;
//     RooRealVar * newPar;
    int i=0;
    cout<<"=============== Original parameters ==============="<<endl;
//    while (ivar) {
    for (int i=0;i<numParameters;++i){
      RooRealVar* ivar =  (RooRealVar* )bkg_ang_params->find(Form("p%02d",i)) ;
      if( !ivar) {
       cout<<"Error: "<<Form("p%02d",i)<<" not found !!!"<<endl;
       exit(1);
      }
      cout <<Form("%s=%f",ivar->GetName(),ivar->getVal())<<endl;
    }
//
    if(!covMatrix) {
     cout<<Form("covMatrix not found: Exit!!!")<<endl;
     exit(0);
    }else{
     cout<<Form("covMatrix Nrows=%d ",covMatrix->GetNrows())<<endl;
    }
    TMatrixDSym *covMatrixSym=(TMatrixDSym *)covMatrix;     
    covMatrixSym->Print("f=  %10.3e  ");
// 
    RooWorkspace wsb("wsb","workspace sideband");
    wsb.import(*covMatrix,Form("covMatrix_bin%d_%d",Q2Bin,RunEra));
    wsb.import(*bkg_ang_pdf,RooFit::RenameAllVariablesExcept(Form("%d",RunEra),"ctK, ctL, phi"));
//     RooArgList* bkg_ang_params = (RooArgList*)bkg_ang_pdf->getParameters(observables);
    wsb.import(*bkg_mass_sb,RooFit::RenameAllVariablesExcept(Form("%d",RunEra),"tagged_mass"));
//    wsb.import(*bkg_mass_sb,RooFit::RenameVariable("slope",Form("slope^{%d}",RunEra)));
    wsb.import(*max_sbl);
    wsb.import(*min_sbr);
    wsb.writeToFile(Form("savesb_%d_b%d.root.Rename%d",RunEra,Q2Bin,i));
    cout<<Form("================= SAVE %d =================",i)<<endl; 
    
//     
    
        
  } 
