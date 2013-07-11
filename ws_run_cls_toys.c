#include "TFile.h"
#include "TString.h"
#include "TROOT.h"
#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooRealVar.h"
#include "RooProfileLL.h"
#include "RooFitResult.h"

#include "RooStats/HypoTestResult.h"
#include "RooStats/HypoTestInverter.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/ToyMCSampler.h"

#include <iostream>
#include <fstream>

using namespace RooFit;
using namespace RooStats;

void ws_run_cls_toys( const char* wsfilename, double testValMin, double testValMax, TString outFile ) {

   TFile* wsfile = TFile::Open( wsfilename ) ;
   if ( wsfile==0x0 ) {
      printf("\n\n\n *** problem opening workspace file %s\n\n", wsfilename ) ;
      return ;
   } else {
      printf("\n\n Opened workspace file %s\n", wsfilename ) ;
   }

   RooWorkspace* w = (RooWorkspace*) wsfile -> Get("ws") ;
   if ( w==0x0 ) {
      printf("\n\n\n *** Did not find workspace in %s\n\n\n", wsfilename ) ;
      return ;
   } else {
      printf(" Found workspace.\n" ) ;
   }

   ModelConfig* mc = (ModelConfig*) w -> obj( "SbModel" ) ;
   if ( mc==0x0 ) {
      printf("\n\n\n *** Did not find SbModel ModelConfig.\n\n") ;
      return ;
   } else {
      printf(" Found SbModel ModelConfig.\n") ;
   }

   RooAbsData* data = w -> data( "hbb_observed_rds" ) ;
   if ( data==0x0 ) {
      printf("\n\n\n *** Did not find RooDataset hbb_observed_rds.\n\n") ;
      return ;
   } else {
      printf(" Found RooDataset hbb_observed_rds.\n") ;
   }

   RooRealVar* myPOI = (RooRealVar*) mc -> GetParametersOfInterest() -> first();
   if ( myPOI==0x0 ) {
      printf("\n\n\n *** didn't find POI in ModelConfig.\n\n\n") ;
      return ;
   } else {
      printf(" Found POI in ModelConfig.\n") ;
   }
   myPOI -> setRange( 0., 100. ) ;

   ModelConfig* bModel = (ModelConfig*) w->obj("BModel");
   ModelConfig* sbModel = (ModelConfig*) w->obj("SbModel");

   ProfileLikelihoodTestStat profll(*sbModel->GetPdf());
   profll.SetOneSided(1);
   TestStatistic* testStat = &profll ;


   //-------------------

   ///AsymptoticCalculator asymptotic_calculator(*data, *bModel, *sbModel);
   ///HypoTestInverter calc_asymptotic_calculator(asymptotic_calculator);

   //-------------------

   HypoTestCalculatorGeneric* hypo_calc = new FrequentistCalculator( *data, *bModel, *sbModel ) ;

   ToyMCSampler* toymcs = (ToyMCSampler*) hypo_calc -> GetTestStatSampler() ;
   toymcs -> SetNEventsPerToy(1) ;
   toymcs -> SetTestStatistic( testStat ) ;

   int nToys = 1000 ;

   ((FrequentistCalculator*) hypo_calc) -> SetToys( nToys, nToys ) ;

   HypoTestInverter calc_toys_calculator( *hypo_calc ) ;

   //-------------------


  calc_toys_calculator.SetConfidenceLevel(0.95);
  calc_toys_calculator.SetVerbose(true);
  calc_toys_calculator.UseCLs(true);



  calc_toys_calculator.SetFixedScan( 5 , testValMin , testValMax) ;
  HypoTestInverterResult* res_toys_calculator = calc_toys_calculator.GetInterval();


//SamplingDistribution * s = res_asymptotic_calculator->GetExpectedPValueDist(0);
//if ( s==0x0 ) {
//   printf("\n\n\n *** Problem getting SamplingDistribution from res_asymptotic_calculator.\n\n") ;
//   return ;
//} else {
//   printf("\n\n Have SamplingDistribution from res_asymptotic_calculator.\n") ;
//}
//const std::vector<double> & values = s->GetSamplingDistribution();

//double maxSigma = res_asymptotic_calculator->fgAsymptoticMaxSigma;
//double dsig = 2* maxSigma/ (values.size() -1) ;

//int  im = (int) TMath::Floor ( ( -1 + maxSigma )/dsig + 0.5);
//int  i0 = (int) TMath::Floor ( ( maxSigma )/dsig + 0.5);
//int  ip = (int) TMath::Floor ( ( +1 + maxSigma )/dsig + 0.5);

//double CLs_expM = values[im];
//double CLs_exp0 = values[i0];
//double CLs_expP = values[ip];



  ofstream outStream ;
  outStream.open(outFile,ios::app) ;

  outStream << "CLs = " << res_toys_calculator->UpperLimit()
    << "   CLs_exp = " << res_toys_calculator->GetExpectedUpperLimit(0)
    << "   CLs_exp(-1s) = " << res_toys_calculator->GetExpectedUpperLimit(-1)
    << "   CLs_exp(+1s) = " << res_toys_calculator->GetExpectedUpperLimit(1) << endl ;


  outStream.close() ;



} // ws_run_cls_toys.



