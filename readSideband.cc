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
#include <RooArgSet.h>
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
char OutSaveFileName[400] =  "savesb_2016_b0.root";
char *ReadOutSaveFileName;

std::string PNGNameReadSB3D="readSideband";
RooDataSet *geneDataNew=0;
RooDataSet *geneMassNew=0;
RooBernsteinSideband *BernSideBand=0;
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

int GenEntries=1000;

int Q2Bin  =0 ;
int RunEra =0 ;

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
    RooRealVar* tagged_mass = new RooRealVar("tagged_mass" , "#mu^{+}#mu^{-}K#pi mass", XMinSign,XMaxSign, "GeV/c^2");
    
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
     cout<<Form("Workspace Found in file : %s. Try to print=>\n",OutSaveFileName)<<endl;
     w->Print();
    }
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
      BernSideBand= (RooBernsteinSideband *)w->pdf(Form("BernSideBand_bin%d_201%d",i,j));
      if ( !BernSideBand) {
//       cout<<Form("BernSideBand_201%d not found\n",i)<<endl;
       if(i==7&&j==8){
//         cout<<Form("BernSideBand_bin%d_201%d not found in file:%s\n",i,j,OutSaveFileName)<<endl;
         cout<<Form("Not found any BernSideBand function in  in file:%s\n",OutSaveFileName)<<endl;
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
//=== dump par    
    int  maxDegree1 =   BernSideBand->maxDegree1Get();
    int  maxDegree2 =   BernSideBand->maxDegree2Get();
    int  maxDegree3 =   BernSideBand->maxDegree3Get();
    int  numParameters = (maxDegree1+1)*(maxDegree2+1)*(maxDegree3+1);

    cout<<Form("maxDegree1=%d\n",maxDegree1); 
    cout<<Form("maxDegree2=%d\n",maxDegree2); 
    cout<<Form("maxDegree3=%d\n",maxDegree3); 
    cout<<Form("numParameters=%d\n",numParameters); 

    RooArgList observables (*ctK, *ctL, *phi);
    RooArgList* bkg_ang_params = (RooArgList*)BernSideBand ->getParameters(observables);
    bool test1 = false;
    bool test2 = false;
    cout<<"Trying to read in the format with 3 digit: "<<Form("p%03d",0)<<endl;
    for (int i=0;i<numParameters;++i){
      RooRealVar* ivar =  (RooRealVar* )bkg_ang_params->find(Form("p%03d",i)) ;
      if( !ivar) {
       cout<<"Error: "<<Form("p%03d",i)<<" not found !!!"<<endl;
       break;
      }else{
       test1 = true;
      } 
      cout <<Form("%s=%f",ivar->GetName(),ivar->getVal())<<endl;
    }
    cout<<"Trying to read in the format with 3 digit and Era: "<<Form("p%03d_%d",0,RunEra)<<endl;
    for (int i=0;i<numParameters;++i){
      RooRealVar* ivar =  (RooRealVar* )bkg_ang_params->find(Form("p%03d_%d",i,RunEra)) ;
      if( !ivar) {
       cout<<"Error: "<<Form("p%03d_%d",i,RunEra)<<" not found !!!"<<endl;
       break;
      }else{
       test2 = true;
      } 
      cout <<Form("%s=%f",ivar->GetName(),ivar->getVal())<<endl;
    }
    if(!test1&&!test2){
     cout<<"Try to dump the params without format"<<endl;
     auto iter  = bkg_ang_params->createIterator();
     RooRealVar* ivar =  (RooRealVar*)iter->Next();
     while (ivar) {
       cout <<Form("%s=%f ",ivar->GetName(),ivar ->getVal())<<endl;
       ivar =  (RooRealVar*)iter->Next();
     }
    } 
//=== dump par    



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
    
    
    RooAbsPdf::GenSpec* genSpec = BernSideBand->prepareMultiGen(RooArgSet(*ctL,*ctK,*phi),NumEvents(GenEntries)) ;
    geneDataNew = BernSideBand->generate(*genSpec) ;
//
    RooAbsPdf::GenSpec* massSpec = bkg_mass_sb->prepareMultiGen(RooArgSet(*tagged_mass ),NumEvents(GenEntries)) ;    
    geneMassNew = bkg_mass_sb->generate(*massSpec) ;									 
    int SBEntries=GenEntries;
														      
    TCanvas* cSideband = new TCanvas("Sideband Gen","Sideband Projections",200,10,900,780);
    cSideband->Divide(2,2);
    RooPlot* xframeFits=ctL->frame(RooFit::Bins(xCosLHBin));
    geneDataNew->plotOn(xframeFits, RooFit::MarkerStyle(kPlus));
    BernSideBand->plotOn( xframeFits,RooFit::LineColor(kRed),Normalization(SBEntries,RooAbsReal::NumEvent),RooFit::ProjWData(RooArgSet(*ctL),*geneDataNew));
    cSideband->cd(1);
    xframeFits->Draw();
    RooPlot* yframeFits=ctK->frame(RooFit::Bins(xCosKHBin));
    geneDataNew->plotOn(yframeFits, RooFit::MarkerStyle(kPlus));
    BernSideBand->plotOn( yframeFits,RooFit::LineColor(kRed),Normalization(SBEntries,RooAbsReal::NumEvent),RooFit::ProjWData(RooArgSet(*ctK),*geneDataNew));
    cSideband->cd(2);
    yframeFits->Draw();
    RooPlot* zframeFits=phi->frame(RooFit::Bins(xPhiHBin));
    geneDataNew->plotOn(zframeFits, RooFit::MarkerStyle(kPlus));
    BernSideBand->plotOn( zframeFits,RooFit::LineColor(kRed),Normalization(SBEntries,RooAbsReal::NumEvent),RooFit::ProjWData(RooArgSet(*phi),*geneDataNew));
    cSideband->cd(3);
    zframeFits->Draw();
    cSideband->cd(4);
    RooPlot* mframeFits=tagged_mass->frame(RooFit::Bins(xMassHBin));
    geneMassNew->plotOn(mframeFits, RooFit::MarkerStyle(kPlus));
    bkg_mass_sb->plotOn(mframeFits,RooFit::LineColor(kRed),Normalization(SBEntries,RooAbsReal::NumEvent),RooFit::ProjWData(RooArgSet(*tagged_mass),*geneMassNew));
    mframeFits->Draw();
    PNGNameReadSB3D=Form("%s-%d-b%d.png",PNGNameReadSB3D.c_str(),RunEra,Q2Bin);

    gSystem->Exec(Form("mv %s %s.tmp",PNGNameReadSB3D.c_str(),PNGNameReadSB3D.c_str()));
    cSideband->Print(PNGNameReadSB3D.c_str());
  } 
