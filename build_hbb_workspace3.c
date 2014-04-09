
#include "TMath.h"
#include "TSystem.h"

#include "RooArgSet.h"
#include "RooConstVar.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooWorkspace.h"
#include "RooPoisson.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooArgSet.h"
#include "RooUniform.h"
#include "RooStats/ModelConfig.h"
#include "TFile.h"
#include "TH1F.h"

#include "RooPosDefCorrGauss.h"
#include "RooAsymAbsProd.h"

#include "getFileValue.c"

#include <fstream>

using namespace RooFit ;
using namespace RooStats ;

char btag_catname[4][10] = { "NT", "2b", "3b", "4b" } ;

const int bins_of_nb(3) ;
const int sigmc_bins_of_nb(4) ;
const int max_bins_of_met(50) ;
int       bins_of_met ;
int       first_met_bin_array_index(0) ;

RooArgSet* globalObservables ;
RooArgSet* allNuisances ;
RooArgSet* allNuisancePdfs ;

int n_shape_systs(0) ;

int syst_type(0) ;

//--- prototypes here.

RooAbsReal* makeLognormalConstraint( const char* NP_name, double NP_val, double NP_err ) ;
RooAbsReal* makeGaussianConstraint( const char* NP_name, double NP_val, double NP_err, bool allowNegative = false ) ;
RooAbsReal* makeCorrelatedLognormalConstraint( const char* NP_name, double NP_val, double NP_err, const char* NP_base_name, bool changeSign=false ) ;
RooAbsReal* makeCorrelatedGaussianConstraint( const char* NP_name, double NP_val, double NP_err, const char* NP_base_name, bool changeSign=false, bool allowNegative = false ) ;

bool setupShapeSyst( const char* infile, const char* systName,
                     int constraintType, // 1=gaussian, 2=...
                     double target_mgl, double target_mlsp,
                     RooWorkspace& workspace
                     ) ;

bool setupBtagSFFracMatrix( const char* infracfile, float sig_mass, const char* postfix, RooWorkspace& workspace ) ;

bool readSignalCounts( const char* susy_counts_file, float sig_mass, RooWorkspace& workspace, const char* zh_counts_file, float chi_to_H_bf ) ;

void makeSystHists( RooWorkspace& workspace, const char* syst_name, const char* syst_var_par, const char* outfile ) ;
void makeSystHistsMultiPar( RooWorkspace& workspace, const char* syst_name, int npars, const char syst_var_pars[][100], const char* outfile ) ;

//===========================================================================================

void build_hbb_workspace3( const char* infile = "outputfiles/input-file.txt",
                           const char* outfile = "outputfiles/ws.root",
                           float sig_mass = 250.,
                           bool use3b = true,
                           bool combine_top_metbins = false,
                           int arg_syst_type = 2, // 1 = Gaussian, 2 = log-normal
                           bool drop_first_met_bin = false,
                           float chi_to_H_bf = 1.0
                           ) {


  //-------------------------------------------------------------------------

  syst_type = arg_syst_type ;

  //-- Create workspace and other RooStats things.

  printf("\n\n Creating workspace.\n\n") ;

  RooWorkspace workspace("ws") ;
  workspace.autoImportClassCode(true) ;

  globalObservables      = new RooArgSet("globalObservables");
  allNuisances           = new RooArgSet("allNuisances");
  allNuisancePdfs        = new RooArgSet("allNuisancePdfs");
  RooArgSet* observedParametersList = new RooArgSet("observables") ;




  //-------------------------------------------------------------------------

  printf("\n\n Reading input file: %s\n\n", infile ) ;

  float fileVal ;
  char pname[1000] ;
  char formula[1000] ;



  sprintf( pname, "bins_of_met" ) ;
  if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
  bins_of_met = TMath::Nint( fileVal ) ;

  //-- save bins_of_met in the workspace for convenience.
  RooRealVar bom( "bins_of_met", "bins_of_met", bins_of_met, 0., 1000. ) ;
  bom.setConstant(kTRUE) ;
  workspace.import(bom) ;

  if ( !drop_first_met_bin ) {
    first_met_bin_array_index = 0 ;
  } else {
    first_met_bin_array_index = 1 ;
  }
  RooRealVar fmbai( "first_met_bin_array_index", "first_met_bin_array_index", first_met_bin_array_index, -1, 2 ) ;
  fmbai.setConstant(kTRUE) ;
  workspace.import(fmbai) ;




  //-- get signal input file and look for requested signal mass.
  char susy_counts_filename[10000] ;
  if ( !getFileStringValue( infile, "signal_counts_file", susy_counts_filename ) ) {
    printf("\n\n *** Can't find input susy counts file: signal_counts_file line of %s.\n\n", infile ) ;
    return ;
  }
  char zh_counts_filename[10000] ;
  if ( !getFileStringValue( infile, "sig_ZH_counts_file", zh_counts_filename ) ) {
    printf("\n\n *** Can't find input susy ZH counts file: sig_ZH_counts_file line of %s.\n\n", infile ) ;
    /// return ; // allow to keep going without this.
    sprintf( zh_counts_filename, "" ) ;
  }


  if ( !readSignalCounts( susy_counts_filename, sig_mass, workspace, zh_counts_filename, chi_to_H_bf ) ) {
    printf("\n\n *** Can't find signal mass of %.0f in %s\n\n", sig_mass, susy_counts_filename ) ;
    return ;
  }






  //-- save bins_of_nb in the workspace for convenience.
  int save_bins_of_nb = 3 ;
  if ( !use3b ) save_bins_of_nb = 2 ;
  RooRealVar bonb( "bins_of_nb", "bins_of_nb", save_bins_of_nb, 0., 1000. ) ;
  bonb.setConstant(kTRUE) ;
  workspace.import(bonb) ;


  RooRealVar* rv_N_msig[bins_of_nb][max_bins_of_met] ; // first index is number of btags, second is met bin.
  RooRealVar* rv_N_msb[bins_of_nb][max_bins_of_met]  ; // first index is number of btags, second is met bin.

  RooAbsReal* rv_Rsigsb_corr[bins_of_nb][max_bins_of_met]  ;

  for ( int nbi=0; nbi<bins_of_nb; nbi++ ) {

    if ( (!use3b) && nbi==1 ) continue ;

    for ( int mbi=first_met_bin_array_index; mbi<bins_of_met; mbi++ ) {

      sprintf( pname, "N_%db_msig_met%d", nbi+2, mbi+1 ) ;
      if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
      rv_N_msig[nbi][mbi] = new RooRealVar( pname, pname, 0., 1.e6 ) ;
      rv_N_msig[nbi][mbi] -> setVal( TMath::Nint(fileVal) ) ;
      rv_N_msig[nbi][mbi] -> setConstant( kTRUE ) ;
      observedParametersList -> add( *rv_N_msig[nbi][mbi] ) ;

      //--- Save input val without Nint so that they can be retrieved later for toy generation means.
      sprintf( pname, "N_%db_msig_met%d_realVal", nbi+2, mbi+1 ) ;
      RooConstVar rcv_msig( pname, pname, fileVal ) ;
      workspace.import( rcv_msig ) ;


      sprintf( pname, "N_%db_msb_met%d", nbi+2, mbi+1 ) ;
      if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
      rv_N_msb[nbi][mbi] = new RooRealVar( pname, pname, 0., 1.e6 ) ;
      rv_N_msb[nbi][mbi] -> setVal( TMath::Nint(fileVal) ) ;
      rv_N_msb[nbi][mbi] -> setConstant( kTRUE ) ;
      observedParametersList -> add( *rv_N_msb[nbi][mbi] ) ;

      //--- Save input val without Nint so that they can be retrieved later for toy generation means.
      sprintf( pname, "N_%db_msb_met%d_realVal", nbi+2, mbi+1 ) ;
      RooConstVar rcv_msb( pname, pname, fileVal ) ;
      workspace.import( rcv_msb ) ;


      if ( (!combine_top_metbins) || mbi==0 ) {

        float corrVal, corrSyst ;
        sprintf( pname, "Rsigsb_syst_%db_met%d", nbi+2, mbi+1 ) ;
        if ( !getFileValue( infile, pname, corrSyst ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
        sprintf( pname, "Rsigsb_corr_%db_met%d", nbi+2, mbi+1 ) ;
        if ( !getFileValue( infile, pname, corrVal  ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
        if ( syst_type == 1 ) {
          rv_Rsigsb_corr[nbi][mbi] = makeGaussianConstraint( pname, corrVal, corrSyst ) ;
        } else if ( syst_type == 2 ) {
          rv_Rsigsb_corr[nbi][mbi] = makeLognormalConstraint( pname, corrVal, corrSyst ) ;
        } else {
          printf("\n\n *** Illegal syst_type %d\n\n", syst_type ) ; return ;
        }

      } else {

        if ( mbi==1 ) {
          float corrVal, corrSyst ;
          sprintf( pname, "Rsigsb_syst_%db_metbins234", nbi+2 ) ;
          if ( !getFileValue( infile, pname, corrSyst ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
          sprintf( pname, "Rsigsb_corr_%db_metbins234", nbi+2 ) ;
          if ( !getFileValue( infile, pname, corrVal  ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
          if ( syst_type == 1 ) {
            rv_Rsigsb_corr[nbi][mbi] = makeGaussianConstraint( pname, corrVal, corrSyst ) ;
          } else if ( syst_type == 2 ) {
            rv_Rsigsb_corr[nbi][mbi] = makeLognormalConstraint( pname, corrVal, corrSyst ) ;
          } else {
            printf("\n\n *** Illegal syst_type %d\n\n", syst_type ) ; return ;
          }
        } else {
          rv_Rsigsb_corr[nbi][mbi] = rv_Rsigsb_corr[nbi][1] ;
        }

      }


    } // mbi.

  } // nbi.


  //-- Get list of shape systs.
  char shape_syst_names[50][1000] ;
  if ( !getFileMultiStringValue( infile, "list_of_shape_systs", n_shape_systs, shape_syst_names ) ) {
    printf("\n\n *** Could not find line starting with list_of_shape_systs in %s\n\n", infile ) ;
    return ;
  }

  for ( int ssi=0; ssi<n_shape_systs; ssi++ ) {
    char shape_syst_file[10000] ;
    char systname[1000] ;
    sprintf( systname, "shape_syst_%s", shape_syst_names[ssi] ) ;
    if ( !getFileStringValue( infile, systname, shape_syst_file ) ) {
      printf("\n\n *** Can't find file for shape syst %s in %s\n\n", shape_syst_names[ssi], infile ) ;
      return ;
    }
    printf("\n\n ======= Reading in shape syst %s from %s\n\n", shape_syst_names[ssi], shape_syst_file ) ;
    setupShapeSyst( shape_syst_file, systname, syst_type, sig_mass, 1., workspace ) ;
    printf("\n back from setupShapeSyst.\n\n" ) ;
    fflush(stdout) ;
  } // ssi.





  //-- Read in the btag SF systematic transfer fraction matrices
  char btag_sf_frac_file[10000] ;
  if ( !getFileStringValue( infile, "btag_SF_frac_matrix_file_SIG", btag_sf_frac_file ) ) {
    printf("\n\n *** Can't find btag SF transfer fraction matrix file name btag_SF_frac_matrix_file_SIG in %s\n\n", infile ) ;
    return ;
  }
  if ( !setupBtagSFFracMatrix( btag_sf_frac_file, sig_mass, "SIG_METsigAll", workspace ) ) {
    printf("\n\n *** Problem reading in btag SF transfer matrix from %s\n\n", btag_sf_frac_file ) ;
    return ;
  }

  if ( !getFileStringValue( infile, "btag_SF_frac_matrix_file_SB", btag_sf_frac_file ) ) {
    printf("\n\n *** Can't find btag SF transfer fraction matrix file name btag_SF_frac_matrix_file_SB in %s\n\n", infile ) ;
    return ;
  }
  if ( !setupBtagSFFracMatrix( btag_sf_frac_file, sig_mass, "SB_METsigAll", workspace ) ) {
    printf("\n\n *** Problem reading in btag SF transfer matrix from %s\n\n", btag_sf_frac_file ) ;
    return ;
  }


  //-- Read in trigger efficiency corrections and uncertainties.

  float trig_eff_corr_val ;
  float trig_eff_corr_err ;

  //-- trig eff, metsig1
  sprintf( pname, "trig_eff_corr_metsig1_val" ) ;
  if ( !getFileValue( infile, pname, trig_eff_corr_val ) ) { printf("\n\n *** error: can't find %s in %s\n\n", pname, infile ) ; return ; }
  sprintf( pname, "trig_eff_corr_metsig1_err" ) ;
  if ( !getFileValue( infile, pname, trig_eff_corr_err ) ) { printf("\n\n *** error: can't find %s in %s\n\n", pname, infile ) ; return ; }

  RooAbsReal* rv_trig_eff_corr_metsig1(0x0) ;
  sprintf( pname, "trig_eff_corr_metsig1" ) ;
  if ( syst_type == 1 ) {
    rv_trig_eff_corr_metsig1 = makeGaussianConstraint( pname, trig_eff_corr_val, trig_eff_corr_err ) ;
  } else if ( syst_type == 2 ) {
    rv_trig_eff_corr_metsig1 = makeLognormalConstraint( pname, trig_eff_corr_val, trig_eff_corr_err ) ;
  } else {
    printf("\n\n *** Illegal syst_type %d\n\n", syst_type ) ; return ;
  }


  //-- trig eff, metsig2
  sprintf( pname, "trig_eff_corr_metsig2_val" ) ;
  if ( !getFileValue( infile, pname, trig_eff_corr_val ) ) { printf("\n\n *** error: can't find %s in %s\n\n", pname, infile ) ; return ; }
  sprintf( pname, "trig_eff_corr_metsig2_err" ) ;
  if ( !getFileValue( infile, pname, trig_eff_corr_err ) ) { printf("\n\n *** error: can't find %s in %s\n\n", pname, infile ) ; return ; }

  RooAbsReal* rv_trig_eff_corr_metsig2(0x0) ;
  sprintf( pname, "trig_eff_corr_metsig2" ) ;
  if ( syst_type == 1 ) {
    rv_trig_eff_corr_metsig2 = makeGaussianConstraint( pname, trig_eff_corr_val, trig_eff_corr_err ) ;
  } else if ( syst_type == 2 ) {
    rv_trig_eff_corr_metsig2 = makeLognormalConstraint( pname, trig_eff_corr_val, trig_eff_corr_err ) ;
  } else {
    printf("\n\n *** Illegal syst_type %d\n\n", syst_type ) ; return ;
  }



  //-- trig eff, metsig34
  sprintf( pname, "trig_eff_corr_metsig34_val" ) ;
  if ( !getFileValue( infile, pname, trig_eff_corr_val ) ) { printf("\n\n *** error: can't find %s in %s\n\n", pname, infile ) ; return ; }
  sprintf( pname, "trig_eff_corr_metsig34_err" ) ;
  if ( !getFileValue( infile, pname, trig_eff_corr_err ) ) { printf("\n\n *** error: can't find %s in %s\n\n", pname, infile ) ; return ; }

  RooAbsReal* rv_trig_eff_corr_metsig34(0x0) ;
  sprintf( pname, "trig_eff_corr_metsig34" ) ;
  if ( syst_type == 1 ) {
    rv_trig_eff_corr_metsig34 = makeGaussianConstraint( pname, trig_eff_corr_val, trig_eff_corr_err ) ;
  } else if ( syst_type == 2 ) {
    rv_trig_eff_corr_metsig34 = makeLognormalConstraint( pname, trig_eff_corr_val, trig_eff_corr_err ) ;
  } else {
    printf("\n\n *** Illegal syst_type %d\n\n", syst_type ) ; return ;
  }




  //-- lumi uncertainty
  RooAbsReal* rv_lumi_syst(0x0) ;
  sprintf( pname, "luminosity_uncertainty" ) ;
  if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** error: can't find %s in %s\n\n", pname, infile ) ; return ; }
  if ( syst_type == 1 ) {
    rv_lumi_syst = makeGaussianConstraint( pname, 1.0, fileVal ) ;
  } else if ( syst_type == 2 ) {
    rv_lumi_syst = makeLognormalConstraint( pname, 1.0, fileVal ) ;
  } else {
    printf("\n\n *** Illegal syst_type %d\n\n", syst_type ) ; return ;
  }


  //-- background sample comp. syst
  RooAbsReal* rv_bg_sample_comp_syst(0x0) ;
  sprintf( pname, "background_sample_comp" ) ;
  if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** error: can't find %s in %s\n\n", pname, infile ) ; return ; }
  if ( syst_type == 1 ) {
    rv_bg_sample_comp_syst = makeGaussianConstraint( pname, 1.0, fileVal ) ;
  } else if ( syst_type == 2 ) {
    rv_bg_sample_comp_syst = makeLognormalConstraint( pname, 1.0, fileVal ) ;
  } else {
    printf("\n\n *** Illegal syst_type %d\n\n", syst_type ) ; return ;
  }







  //-- Finished reading input from file.

  //-------------------------------------------------------------------------

  printf("\n\n Creating and importing dataset into workspace.\n\n") ;

  RooDataSet* dsObserved = new RooDataSet("hbb_observed_rds", "hbb observed data values", *observedParametersList ) ;
  dsObserved -> add( *observedParametersList ) ;
  workspace.import( *dsObserved ) ;

  //-------------------------------------------------------------------------

  //-- Define all floats.

  printf("\n\n Defining all unconstrained floats (Ratios, signal strength).\n\n") ;

  double R_msigmsb_initialval(0.15) ;

  RooRealVar* rv_R_msigmsb[50] ;

  for ( int mbi=first_met_bin_array_index; mbi<bins_of_met; mbi++ ) {

    if ( (!combine_top_metbins) || mbi==0 ) {

      sprintf( pname, "R_msigmsb_met%d", mbi+1 ) ;
      printf( "  %s\n", pname ) ;
      rv_R_msigmsb[mbi] = new RooRealVar( pname, pname, R_msigmsb_initialval, 0., 3. ) ;
      rv_R_msigmsb[mbi] -> setConstant( kFALSE ) ;
      rv_R_msigmsb[mbi] -> Print() ;

    } else {

      if ( mbi==1 ) {
        sprintf( pname, "R_msigmsb_metbins234" ) ;
        printf( "  %s\n", pname ) ;
        rv_R_msigmsb[mbi] = new RooRealVar( pname, pname, R_msigmsb_initialval, 0., 3. ) ;
        rv_R_msigmsb[mbi] -> setConstant( kFALSE ) ;
        rv_R_msigmsb[mbi] -> Print() ;
      } else {
        rv_R_msigmsb[mbi] = rv_R_msigmsb[1] ;
      }

    }

  } // mbi.

  printf("\n") ;

  sprintf( pname, "sig_strength" ) ;
  RooRealVar* rv_sig_strength = new RooRealVar( pname, pname, 1.0, 0., 10. ) ;
  rv_sig_strength -> setConstant(kFALSE) ;
  rv_sig_strength -> Print() ;
  printf("  %s\n\n", pname ) ;

  //-------------------------------------------------------------------------

  //-- Define mu parameters for background.

  printf("\n\n Defining mu parameters for background.\n\n") ;

  RooAbsReal* rv_mu_bg_msig[bins_of_nb][max_bins_of_met] ;  // first index is number of btags, second is met bin.
  RooAbsReal* rv_mu_bg_msb[bins_of_nb][max_bins_of_met]  ;  // first index is number of btags, second is met bin.

  for ( int nbi=0; nbi<bins_of_nb; nbi++ ) {

    if ( (!use3b) && nbi==1 ) continue ;

    for ( int mbi=first_met_bin_array_index; mbi<bins_of_met; mbi++ ) {

      sprintf( pname, "mu_bg_%db_msb_met%d", nbi+2, mbi+1 ) ;
      printf( "  %s\n", pname ) ;
      rv_mu_bg_msb[nbi][mbi] = new RooRealVar( pname, pname, rv_N_msb[nbi][mbi] -> getVal(), 0., 1.e6 ) ;
      rv_mu_bg_msb[nbi][mbi] -> Print() ;


      sprintf( formula, "@0 * @1 * @2 * @3" ) ;
      sprintf( pname, "mu_bg_%db_msig_met%d", nbi+2, mbi+1 ) ;
      printf( "  %s\n", pname ) ;
      rv_mu_bg_msig[nbi][mbi] = new RooFormulaVar( pname, formula, RooArgSet( *rv_bg_sample_comp_syst, *rv_Rsigsb_corr[nbi][mbi], *rv_R_msigmsb[mbi], *rv_mu_bg_msb[nbi][mbi] ) ) ;
      rv_mu_bg_msig[nbi][mbi] -> Print() ;


    } // mbi.

  } // nbi.





  //-- Define mu parameters for signal.

  printf("\n\n Defining mu parameters for signal.\n\n") ;

  RooAbsReal* rv_mu_sig_msig_nobtagsyst[sigmc_bins_of_nb][max_bins_of_met] ; // first index is number of btags, second is met bin.
  RooAbsReal* rv_mu_sig_msb_nobtagsyst[sigmc_bins_of_nb][max_bins_of_met]  ; // first index is number of btags, second is met bin.


  for ( int smcnbi=0; smcnbi<sigmc_bins_of_nb; smcnbi++ ) {  //--- this is 0 = no-tag, 1 = 2b, 2 = 3b, 3 = 4b

    for ( int mbi=first_met_bin_array_index; mbi<bins_of_met; mbi++ ) {

      //-- set up combination of all signal shape systematics.
      char syst_prod_eqn[1000] ;
      sprintf( syst_prod_eqn, "@0" ) ;
      for ( int ssi=1; ssi<n_shape_systs; ssi++ ) {
        char tmpstr[1000] ;
        sprintf( tmpstr, "%s * @%d", syst_prod_eqn, ssi ) ;
        sprintf( syst_prod_eqn, "%s", tmpstr ) ;
      } // ssi.

      RooArgSet shapeSystProdSet_msig ;
      for ( int ssi=0; ssi<n_shape_systs; ssi++ ) {
        sprintf( pname, "shape_syst_%s_msig_met%d_%s", shape_syst_names[ssi], mbi+1, btag_catname[smcnbi] ) ;
        RooAbsReal* rar_sf = (RooAbsReal*) workspace.obj( pname ) ;
        if ( rar_sf == 0x0 ) { printf("\n\n *** Missing %s shape syst for met %d, nb %s,  (%s)\n\n", shape_syst_names[ssi], mbi+1, btag_catname[smcnbi], pname ) ; return ; }
        shapeSystProdSet_msig.add( *rar_sf ) ;
      } // ssi.
      sprintf( pname, "shape_syst_prod_msig_met%d_%s", mbi+1, btag_catname[smcnbi] ) ;
      RooFormulaVar* rfv_shape_syst_prod_msig = new RooFormulaVar( pname, syst_prod_eqn, shapeSystProdSet_msig ) ;

      RooArgSet shapeSystProdSet_msb ;
      for ( int ssi=0; ssi<n_shape_systs; ssi++ ) {
        sprintf( pname, "shape_syst_%s_msb_met%d_%s", shape_syst_names[ssi], mbi+1, btag_catname[smcnbi] ) ;
        RooAbsReal* rar_sf = (RooAbsReal*) workspace.obj( pname ) ;
        if ( rar_sf == 0x0 ) { printf("\n\n *** Missing %s shape syst for met %d, nb %s,  (%s)\n\n", shape_syst_names[ssi], mbi+1, btag_catname[smcnbi], pname ) ; return ; }
        shapeSystProdSet_msb.add( *rar_sf ) ;
      } // ssi.
      sprintf( pname, "shape_syst_prod_msb_met%d_%s", mbi+1, btag_catname[smcnbi] ) ;
      RooFormulaVar* rfv_shape_syst_prod_msb = new RooFormulaVar( pname, syst_prod_eqn, shapeSystProdSet_msb ) ;



      //-- retrieve signals pars from workspace, created in readSignalCounts

      sprintf( pname, "smc_%s_msig_met%d", btag_catname[smcnbi], mbi+1 ) ;
      RooRealVar* smc_msig = (RooRealVar*) workspace.obj( pname ) ;
      if ( smc_msig == 0x0 ) { printf("\n\n *** missing %s in workspace.\n\n", pname ) ; return ; }

      sprintf( pname, "smc_%s_msb_met%d", btag_catname[smcnbi], mbi+1 ) ;
      RooRealVar* smc_msb = (RooRealVar*) workspace.obj( pname ) ;
      if ( smc_msb == 0x0 ) { printf("\n\n *** missing %s in workspace.\n\n", pname ) ; return ; }

      sprintf( pname, "syst_sig_eff_mc_stats_msig_met%d_%s", mbi+1, btag_catname[smcnbi] ) ;
      RooRealVar* smc_msig_mcstat_syst = (RooRealVar*) workspace.obj( pname ) ;
      if ( smc_msig_mcstat_syst == 0x0 ) { printf("\n\n *** missing %s in workspace.\n\n", pname ) ; return ; }

      sprintf( pname, "syst_sig_eff_mc_stats_msb_met%d_%s", mbi+1, btag_catname[smcnbi] ) ;
      RooRealVar* smc_msb_mcstat_syst = (RooRealVar*) workspace.obj( pname ) ;
      if ( smc_msb_mcstat_syst == 0x0 ) { printf("\n\n *** missing %s in workspace.\n\n", pname ) ; return ; }

      RooAbsReal* rv_trig_eff_corr(0x0) ;
      if ( mbi == 0 ) {
        rv_trig_eff_corr = rv_trig_eff_corr_metsig1 ;
      } else if ( mbi == 1 ) {
        rv_trig_eff_corr = rv_trig_eff_corr_metsig2 ;
      } else {
        rv_trig_eff_corr = rv_trig_eff_corr_metsig34 ;
      }

      sprintf( formula, "@0 * @1 * @2 * @3 * @4 * @5" ) ;
      sprintf( pname, "mu_sig_%s_msig_met%d_nobtagsyst", btag_catname[smcnbi], mbi+1 ) ;
      printf( "  %s, smcnbi=%d, mbi=%d \n", pname, smcnbi, mbi ) ;
      rv_mu_sig_msig_nobtagsyst[smcnbi][mbi] = new RooFormulaVar( pname, formula, RooArgSet( *rv_lumi_syst, *rv_trig_eff_corr, *rv_sig_strength, *rfv_shape_syst_prod_msig, *smc_msig_mcstat_syst, *smc_msig ) ) ;
      rv_mu_sig_msig_nobtagsyst[smcnbi][mbi] -> Print() ;

      sprintf( formula, "@0 * @1 * @2 * @3 * @4 * @5" ) ;
      sprintf( pname, "mu_sig_%s_msb_met%d_nobtagsyst", btag_catname[smcnbi], mbi+1 ) ;
      printf( "  %s, smcnbi=%d, mbi=%d\n", pname, smcnbi, mbi ) ;
      rv_mu_sig_msb_nobtagsyst[smcnbi][mbi] = new RooFormulaVar( pname, formula, RooArgSet( *rv_lumi_syst, *rv_trig_eff_corr, *rv_sig_strength, *rfv_shape_syst_prod_msb, *smc_msb_mcstat_syst, *smc_msb ) ) ;
      rv_mu_sig_msb_nobtagsyst[smcnbi][mbi] -> Print() ;


    } // mbi.

  } // smcnbi.



  RooAbsReal* rv_mu_sig_msig[sigmc_bins_of_nb][max_bins_of_met] ; // first index is number of btags, second is met bin.
  RooAbsReal* rv_mu_sig_msb[sigmc_bins_of_nb][max_bins_of_met]  ; // first index is number of btags, second is met bin.

  for ( int mbi=first_met_bin_array_index; mbi<bins_of_met; mbi++ ) {

    for ( int to_smcnbi=0; to_smcnbi<sigmc_bins_of_nb; to_smcnbi++ ) {  //--- this is 0 = no-tag, 1 = 2b, 2 = 3b, 3 = 4b

      RooArgSet formulaArgs_msig ;
      RooArgSet formulaArgs_msb ;

      for ( int from_smcnbi=0; from_smcnbi<sigmc_bins_of_nb; from_smcnbi++ ) {  //--- this is 0 = no-tag, 1 = 2b, 2 = 3b, 3 = 4b

        sprintf( pname, "frac_from_%s_to_%s_%s", btag_catname[from_smcnbi], btag_catname[to_smcnbi], "SIG_METsigAll" ) ;
        RooAbsReal* frac_sig = (RooAbsReal*) workspace.obj( pname ) ;
        if ( frac_sig == 0x0 ) { printf("\n\n *** missing fraction parameter %s\n\n", pname ) ; return ; }

        sprintf( pname, "frac_from_%s_to_%s_%s", btag_catname[from_smcnbi], btag_catname[to_smcnbi], "SB_METsigAll" ) ;
        RooAbsReal* frac_sb  = (RooAbsReal*) workspace.obj( pname ) ;
        if ( frac_sb  == 0x0 ) { printf("\n\n *** missing fraction parameter %s\n\n", pname ) ; return ; }

        formulaArgs_msig.add( *frac_sig ) ;
        formulaArgs_msig.add( *rv_mu_sig_msig_nobtagsyst[from_smcnbi][mbi] ) ;

        formulaArgs_msb.add(  *frac_sb ) ;
        formulaArgs_msb.add(  *rv_mu_sig_msb_nobtagsyst[from_smcnbi][mbi] ) ;

      } // from_smcnbi.

      sprintf( formula, "@0 * @1 + @2 * @3 + @4 * @5 + @6 * @7" ) ;
      sprintf( pname, "mu_sig_%s_msig_met%d", btag_catname[to_smcnbi], mbi+1 ) ;
      rv_mu_sig_msig[to_smcnbi][mbi] = new RooFormulaVar( pname, formula, formulaArgs_msig ) ;
      rv_mu_sig_msig[to_smcnbi][mbi] -> Print() ;
      if ( to_smcnbi == 0 ) {
        RooAbsReal* rar = (RooAbsReal*) workspace.obj( rv_mu_sig_msig[to_smcnbi][mbi] -> GetName() ) ;
        if ( rar == 0x0 ) {
          printf("\n importing %s into workspace.\n", rv_mu_sig_msig[to_smcnbi][mbi] -> GetName() ) ;
          workspace.import( *rv_mu_sig_msig[to_smcnbi][mbi], RecycleConflictNodes() ) ; // this isn't used later so import it here.
        }
      }


      RooRealVar* rrv_btag_SF_base_par = (RooRealVar*) workspace.obj( "prim_btag_SF_syst" ) ;
      if ( rrv_btag_SF_base_par == 0x0 ) { printf("\n\n *** can't find prim_btag_SF_syst\n\n" ) ; return ; }

      {
        rrv_btag_SF_base_par -> setVal( 0. ) ;
        double val_nom = rv_mu_sig_msig[to_smcnbi][mbi] -> getVal() ;

        rrv_btag_SF_base_par -> setVal( 1. ) ;
        double val_p1s = rv_mu_sig_msig[to_smcnbi][mbi] -> getVal() ;

        rrv_btag_SF_base_par -> setVal(-1. ) ;
        double val_m1s = rv_mu_sig_msig[to_smcnbi][mbi] -> getVal() ;

        rrv_btag_SF_base_par -> setVal( 0. ) ;

        printf("     Test of %s (-1s, nom, +1s) : %.2f, %.2f, %.2f\n\n", pname, val_m1s, val_nom, val_p1s ) ;
      }


      sprintf( formula, "@0 * @1 + @2 * @3 + @4 * @5 + @6 * @7" ) ;
      sprintf( pname, "mu_sig_%s_msb_met%d", btag_catname[to_smcnbi], mbi+1 ) ;
      rv_mu_sig_msb[to_smcnbi][mbi] = new RooFormulaVar( pname, formula, formulaArgs_msb ) ;
      rv_mu_sig_msb[to_smcnbi][mbi] -> Print() ;
      if ( to_smcnbi == 0 ) {
        RooAbsReal* rar = (RooAbsReal*) workspace.obj( rv_mu_sig_msb[to_smcnbi][mbi] -> GetName() ) ;
        if ( rar == 0x0 ) {
          printf("\n importing %s into workspace.\n", rv_mu_sig_msb[to_smcnbi][mbi] -> GetName() ) ;
          workspace.import( *rv_mu_sig_msb[to_smcnbi][mbi], RecycleConflictNodes() ) ; // this isn't used later so import it here.
        }
      }

      {
        rrv_btag_SF_base_par -> setVal( 0. ) ;
        double val_nom = rv_mu_sig_msb[to_smcnbi][mbi] -> getVal() ;

        rrv_btag_SF_base_par -> setVal( 1. ) ;
        double val_p1s = rv_mu_sig_msb[to_smcnbi][mbi] -> getVal() ;

        rrv_btag_SF_base_par -> setVal(-1. ) ;
        double val_m1s = rv_mu_sig_msb[to_smcnbi][mbi] -> getVal() ;

        rrv_btag_SF_base_par -> setVal( 0. ) ;

        printf("     Test of %s (-1s, nom, +1s) : %.2f, %.2f, %.2f\n\n", pname, val_m1s, val_nom, val_p1s ) ;
      }

    } // to_smcnbi.

  } // mbi.

  //-- Finished defining mu parameters.

  //-------------------------------------------------------------------------

  //-- Defining small n's

  printf("\n\n Defining small n's.\n\n") ;

  RooAbsReal* rv_n_msig[bins_of_nb][max_bins_of_met] ;  // first index is number of btags, second is met bin.
  RooAbsReal* rv_n_msb[bins_of_nb][max_bins_of_met]  ;  // first index is number of btags, second is met bin.

  for ( int nbi=0; nbi<bins_of_nb; nbi++ ) {

    if ( (!use3b) && nbi==1 ) continue ;

    int smcnbi = nbi+1 ;

    for ( int mbi=first_met_bin_array_index; mbi<bins_of_met; mbi++ ) {

      sprintf( formula, "@0 + @1" ) ;

      sprintf( pname, "n_%db_msig_met%d", nbi+2, mbi+1 ) ;
      printf( "  %s\n", pname ) ;
      rv_n_msig[nbi][mbi] = new RooFormulaVar( pname, formula, RooArgSet( *rv_mu_sig_msig[smcnbi][mbi], *rv_mu_bg_msig[nbi][mbi] ) ) ;
      rv_n_msig[nbi][mbi] -> Print() ;
      //  workspace.import( *rv_n_msig[nbi][mbi] ) ;

      sprintf( pname, "n_%db_msb_met%d", nbi+2, mbi+1 ) ;
      printf( "  %s\n", pname ) ;
      rv_n_msb[nbi][mbi] = new RooFormulaVar( pname, formula, RooArgSet( *rv_mu_sig_msb[smcnbi][mbi], *rv_mu_bg_msb[nbi][mbi] ) ) ;
      rv_n_msb[nbi][mbi] -> Print() ;
      //  workspace.import( *rv_n_msb[nbi][mbi] ) ;

    } // mbi.

  } // nbi.

    //-------------------------------------------------------------------------

  //-- Define the Poisson pdfs for the observables.

  printf("\n\n Defining Poisson pdfs for the observables.\n\n") ;

  RooAbsReal* rv_pdf_msig[bins_of_nb][max_bins_of_met] ;  // first index is number of btags, second is met bin.
  RooAbsReal* rv_pdf_msb[bins_of_nb][max_bins_of_met]  ;  // first index is number of btags, second is met bin.

  RooArgSet pdflist ;

  for ( int nbi=0; nbi<bins_of_nb; nbi++ ) {

    if ( (!use3b) && nbi==1 ) continue ;

    for ( int mbi=first_met_bin_array_index; mbi<bins_of_met; mbi++ ) {

      sprintf( pname, "pdf_%db_msig_met%d", nbi+2, mbi+1 ) ;
      printf( "  %s\n", pname ) ;
      rv_pdf_msig[nbi][mbi] = new RooPoisson( pname, pname, *rv_N_msig[nbi][mbi], *rv_n_msig[nbi][mbi] ) ;
      rv_pdf_msig[nbi][mbi] -> Print() ;

      pdflist.add( *rv_pdf_msig[nbi][mbi] ) ;

      sprintf( pname, "pdf_%db_msb_met%d", nbi+2, mbi+1 ) ;
      printf( "  %s\n", pname ) ;
      rv_pdf_msb[nbi][mbi] = new RooPoisson( pname, pname, *rv_N_msb[nbi][mbi], *rv_n_msb[nbi][mbi] ) ;
      rv_pdf_msb[nbi][mbi] -> Print() ;

      pdflist.add( *rv_pdf_msb[nbi][mbi] ) ;

    } // mbi.

  } // nbi.

    //-------------------------------------------------------------------------

  //-- Build the likelihood.

  printf("\n\n Building the likelihood.\n\n") ;

  pdflist.add( *allNuisancePdfs ) ;

  pdflist.Print() ;
  printf("\n") ;

  RooProdPdf* likelihood = new RooProdPdf( "likelihood", "hbb likelihood", pdflist ) ;
  likelihood->Print() ;


  //-------------------------------------------------------------------------





  //  printf("\n\n Running a test fit.\n\n") ;

  //  dsObserved -> Print() ;
  //  dsObserved -> printMultiline(cout, 1, kTRUE, "") ;

  //  printf("\n\n =============================================\n\n") ;
  //  likelihood -> fitTo( *dsObserved, PrintLevel(3), Hesse(0), Minos(0) ) ;
  //  printf("\n\n =============================================\n\n") ;







  //-- Set up RooStats models.

  printf("\n\n Setting up S+B model.\n\n") ;

  RooArgSet poi( *rv_sig_strength, "poi" ) ;
  RooUniform signal_prior( "signal_prior", "signal_prior", *rv_sig_strength ) ;

  ModelConfig sbModel ("SbModel");
  sbModel.SetWorkspace( workspace ) ;
  sbModel.SetPdf( *likelihood ) ;
  sbModel.SetParametersOfInterest( poi );
  sbModel.SetPriorPdf(signal_prior);
  sbModel.SetObservables( *observedParametersList );
  sbModel.SetNuisanceParameters( *allNuisances );
  sbModel.SetGlobalObservables( *globalObservables );

  workspace.Print() ;

  printf("\n\n Doing fit for S+B model.\n" ) ; fflush(stdout) ;

  RooAbsReal* pNll = sbModel.GetPdf()->createNLL(*dsObserved);
  RooAbsReal* pProfile = pNll->createProfile(RooArgSet());
  pProfile->getVal();
  RooArgSet* pPoiAndNuisance = new RooArgSet();
  pPoiAndNuisance->add(*sbModel.GetParametersOfInterest());
  if(sbModel.GetNuisanceParameters()) pPoiAndNuisance->add(*sbModel.GetNuisanceParameters());
  printf("\n\n Will save these parameter points that correspond to the fit to data.\n\n") ; fflush(stdout) ;
  pPoiAndNuisance->Print("v");
  sbModel.SetSnapshot(*pPoiAndNuisance);
  workspace.import (sbModel);

  delete pProfile ;
  delete pNll ;
  delete pPoiAndNuisance ;

  printf("\n\n Setting up BG-only model.\n\n") ;

  ModelConfig bModel (*(RooStats::ModelConfig *)workspace.obj("SbModel"));
  bModel.SetName("BModel");
  bModel.SetWorkspace(workspace);

  printf("\n\n Doing fit for BG-only model.\n" ) ; fflush(stdout) ;
  pNll = bModel.GetPdf()->createNLL(*dsObserved);
  pProfile = pNll->createProfile(*bModel.GetParametersOfInterest());
  ((RooRealVar *)(bModel.GetParametersOfInterest()->first()))->setVal(0.);
  pProfile->getVal();
  pPoiAndNuisance = new RooArgSet();
  pPoiAndNuisance->add(*bModel.GetParametersOfInterest());
  if(bModel.GetNuisanceParameters()) pPoiAndNuisance->add(*bModel.GetNuisanceParameters());
  printf("\n\n Should use these parameter points to generate pseudo data for bkg only.\n\n") ; fflush(stdout) ;
  pPoiAndNuisance->Print("v");
  bModel.SetSnapshot(*pPoiAndNuisance);
  workspace.import (bModel);

  delete pProfile ;
  delete pNll ;
  delete pPoiAndNuisance ;

  workspace.Print() ;

  printf("\n\n Saving workspace in : %s\n\n", outfile ) ;

  gSystem->Exec(" mkdir -p outputfiles " ) ;

  workspace.writeToFile( outfile ) ;




  makeSystHists( workspace, "btagSF", "prim_btag_SF_syst", outfile ) ;
  makeSystHists( workspace, "lumi", "prim_luminosity_uncertainty", outfile ) ;

  for ( int ssi=0; ssi<n_shape_systs; ssi++ ) {
    char spname[100] ;
    sprintf( spname, "prim_shape_syst_%s", shape_syst_names[ssi] ) ;
    makeSystHists( workspace, shape_syst_names[ssi], spname, outfile ) ;
  } // ssi.


  int npars(0) ;
  char syst_var_pars[100][100] ;

  sprintf( syst_var_pars[npars++], "prim_trig_eff_corr_metsig1" ) ;
  sprintf( syst_var_pars[npars++], "prim_trig_eff_corr_metsig2" ) ;
  sprintf( syst_var_pars[npars++], "prim_trig_eff_corr_metsig34" ) ;
  makeSystHistsMultiPar( workspace, "trig", npars, syst_var_pars, outfile ) ;

  npars = 0 ;
  for ( int mbi=first_met_bin_array_index; mbi<bins_of_met; mbi++ ) {
    for ( int smcnbi=0; smcnbi<sigmc_bins_of_nb; smcnbi++ ) {  //--- this is 0 = no-tag, 1 = 2b, 2 = 3b, 3 = 4b
      sprintf( syst_var_pars[npars++], "prim_syst_sig_eff_mc_stats_msig_met%d_%s", mbi+1, btag_catname[smcnbi] ) ;
      sprintf( syst_var_pars[npars++], "prim_syst_sig_eff_mc_stats_msb_met%d_%s", mbi+1, btag_catname[smcnbi] ) ;
    } // smcnbi
  } // mbi
  makeSystHistsMultiPar( workspace, "mcstat", npars, syst_var_pars, outfile ) ;



} // build_hbb_workspace3.







  //==============================================================================================

RooAbsReal* makeGaussianConstraint( const char* NP_name, double NP_val, double NP_err, bool allowNegative ) {

  if ( NP_err <= 0. ) {
    printf(" makeGaussianConstraint:  Uncertainty is zero.  Will return constant scale factor of %g for %s.  Input val = %g, err = %g.\n", NP_val, NP_name, NP_val, NP_err ) ;
    return new RooConstVar( NP_name, NP_name, NP_val ) ;
  }


  double max = NP_val + 6.*NP_err ;
  double min = NP_val - 6.*NP_err ;

  if ( min < 0. && !allowNegative ) { min = 1e-5 ; }

  RooRealVar* np_rrv = new RooRealVar( NP_name, NP_name, min, max ) ;
  np_rrv -> setVal( NP_val ) ;
  np_rrv -> setConstant( kFALSE ) ;

  //-- create const variables for mean and sigma so that they can be saved and accessed from workspace later.

  char vname[1000] ;
  sprintf( vname, "mean_%s", NP_name ) ;
  RooRealVar* g_mean = new RooRealVar( vname, vname, NP_val, -1000., 1000. ) ;
  g_mean->setConstant(kTRUE);
  sprintf( vname, "sigma_%s", NP_name ) ;
  RooConstVar* g_sigma = new RooConstVar( vname, vname, NP_err ) ;

  char pdfname[1000] ;
  sprintf( pdfname, "pdf_%s", NP_name ) ;
  RooGaussian* np_pdf = new RooGaussian( pdfname, pdfname, *np_rrv, *g_mean, *g_sigma ) ;

  allNuisances -> add( *np_rrv ) ;
  allNuisancePdfs -> add( *np_pdf ) ;
  globalObservables -> add( *g_mean ) ;

  printf("  makeGaussianConstraint : created nuisance parameter %s : val = %g\n", NP_name, np_rrv -> getVal() ) ;

  return np_rrv ;


} // makeGaussianConstraint.








//==============================================================================================================

RooAbsReal* makeLognormalConstraint( const char* NP_name, double NP_val, double NP_err ) {

  if ( NP_err <= 0. ) {
    printf(" makeLognormalConstraint:  Uncertainty is zero.  Will return constant scale factor of %g for %s.  Input val = %g, err = %g.\n", NP_val, NP_name, NP_val, NP_err ) ;
    return new RooConstVar( NP_name, NP_name, NP_val ) ;
  }

  char pname[1000] ;
  sprintf( pname, "prim_%s", NP_name ) ;

  printf(" makeLognormalConstraint : creating primary log-normal variable %s\n", pname ) ;
  RooRealVar* np_prim_rrv = new RooRealVar( pname, pname, 0., -6., 6. ) ;
  np_prim_rrv -> setVal( 0. ) ;
  np_prim_rrv -> setConstant( kFALSE ) ;

  sprintf( pname, "prim_mean_%s", NP_name ) ;
  RooRealVar* np_prim_mean = new RooRealVar( pname, pname, 0., -6., 6. ) ;
  np_prim_mean->setConstant(kTRUE) ;

  sprintf( pname, "prim_sigma_%s", NP_name ) ;
  RooConstVar* np_prim_sigma = new RooConstVar( pname, pname, 1. ) ;


  char pdfname[1000] ;
  sprintf( pdfname, "pdf_prim_%s", NP_name ) ;
  RooGaussian* np_prim_pdf = new RooGaussian( pdfname, pdfname, *np_prim_rrv, *np_prim_mean, *np_prim_sigma ) ;

  allNuisances -> add( *np_prim_rrv ) ;
  allNuisancePdfs -> add( *np_prim_pdf ) ;
  globalObservables -> add( *np_prim_mean ) ;


  //-- create const variables for mean and sigma so that they can be saved and accessed from workspace later.

  char vname[1000] ;
  sprintf( vname, "mean_%s", NP_name ) ;
  RooConstVar* g_mean  = new RooConstVar( vname, vname, NP_val ) ;
  sprintf( vname, "sigma_%s", NP_name ) ;
  RooConstVar* g_sigma = new RooConstVar( vname, vname, NP_err ) ;

  //-- compute the log-normal-distributed parameter from the primary parameter.

  //--- This is the new way.  RMS of lognormal is much closer to sigma when sigma is
  //    large, doing it this way.  When sigma/mean is small, they are about the same.
  //    That is, exp(sigma/mean) is close to (sigma/mean + 1).  This one is better when
  //    sigma/mean is not small.  The high-side tail is not as strong.
  //
  RooFormulaVar* np_rfv = new RooFormulaVar( NP_name, "@0 * pow( ( @1/@0 + 1. ), @2)",
                                             RooArgSet( *g_mean, *g_sigma, *np_prim_rrv ) ) ;
  //------------------------------------------------------------------------------------------


  printf("  makeLognormalConstraint : created log-normal nuisance parameter %s : val = %g\n", NP_name, np_rfv -> getVal() ) ;

  return np_rfv ;


} // makeLognormalConstraint.







  //==============================================================================================

RooAbsReal* makeCorrelatedGaussianConstraint(
                                             const char* NP_name, double NP_val, double NP_err, const char* NP_base_name, bool changeSign, bool allowNegative ) {

  if ( NP_err <= 0. ) {
    printf("  makeCorrelatedGaussianConstraint: Uncertainty is zero.  Will return constant scale factor of %g for %s.  Input val = %g, err = %g.\n", NP_val, NP_name, NP_val, NP_err ) ;
    return new RooConstVar( NP_name, NP_name, NP_val ) ;
  }

  RooRealVar* rrv_np_base_par = (RooRealVar*) allNuisances -> find( NP_base_name ) ;

  if ( rrv_np_base_par == 0x0 ) {

    printf("\n\n makeCorrelatedGaussianConstraint : creating base nuisance parameter - %s\n\n", NP_base_name ) ;
    rrv_np_base_par = new RooRealVar( NP_base_name, NP_base_name, -6.0, 6.0 ) ;
    rrv_np_base_par -> setVal( 0. ) ;
    rrv_np_base_par -> setConstant( kFALSE ) ;
    allNuisances -> add( *rrv_np_base_par ) ;

    char vname[1000] ;
    sprintf( vname, "mean_%s", NP_base_name ) ;
    RooRealVar* g_mean = new RooRealVar( vname, vname, 0.0,-1000.,1000. ) ;
    g_mean->setConstant(kTRUE);
    sprintf( vname, "sigma_%s", NP_base_name ) ;
    RooConstVar* g_sigma = new RooConstVar( vname, vname, 1.0 ) ;

    char pdfname[100] ;
    sprintf( pdfname, "pdf_%s", NP_base_name ) ;
    printf("\n\n makeCorrelatedGaussianConstraint : creating base nuisance parameter pdf - %s\n\n", pdfname ) ;
    RooGaussian* base_np_pdf = new RooGaussian( pdfname, pdfname, *rrv_np_base_par, *g_mean, *g_sigma ) ;
    allNuisancePdfs -> add( *base_np_pdf ) ;
    globalObservables -> add( *g_mean ) ;

  }

  //-- create const variables for mean and sigma so that they can be saved and accessed from workspace later.

  char vname[1000] ;
  sprintf( vname, "mean_%s", NP_name ) ;
  RooConstVar* g_mean = new RooConstVar( vname, vname, NP_val ) ;
  sprintf( vname, "sigma_%s", NP_name ) ;
  RooConstVar* g_sigma = new RooConstVar( vname, vname, NP_err ) ;

  RooAbsReal* rar(0x0) ;

  if ( allowNegative ) {

    char formula[1000] ;

    if ( !changeSign ) {
      sprintf( formula, "@0+@1*@2" ) ;
    } else {
      sprintf( formula, "@0-@1*@2" ) ;
    }

    rar = new RooFormulaVar( NP_name, formula, RooArgSet( *g_mean, *g_sigma, *rrv_np_base_par ) ) ;

    printf(" makeCorrelatedGaussianConstraint : creating correlated gaussian NP with formula : %s,  %s, val = %g\n", formula, NP_name, rar->getVal() ) ;

  } else {

    rar = new RooPosDefCorrGauss( NP_name, NP_name, *g_mean, *g_sigma, *rrv_np_base_par, changeSign ) ;

    printf(" makeCorrelatedGaussianConstraint : creating pos-def correlated gaussian NP  :  %s, val = %g, err = %g\n", NP_name, rar->getVal(), NP_err ) ;

  }




  return rar ;

} // makeCorrelatedGaussianConstraint.








//==============================================================================================================

RooAbsReal* makeCorrelatedLognormalConstraint(
                                              const char* NP_name, double NP_val, double NP_err, const char* NP_base_name, bool changeSign ) {


  if ( NP_err <= 0. ) {
    printf("  makeCorrelatedLognormalConstraint: Uncertainty is zero.  Will return constant scale factor of %g for %s.  Input val = %g, err = %g.\n", NP_val, NP_name, NP_val, NP_err ) ;
    return new RooConstVar( NP_name, NP_name, NP_val ) ;
  }

  char prim_name[1000] ;
  sprintf( prim_name, "prim_%s", NP_base_name ) ;
  RooRealVar* rrv_np_base_par = (RooRealVar*) allNuisances -> find( prim_name ) ;

  if ( rrv_np_base_par == 0x0 ) {

    printf("\n\n makeCorrelatedLognormalConstraint : creating base nuisance parameter - %s\n\n", prim_name ) ;
    rrv_np_base_par = new RooRealVar( prim_name, prim_name, -6.0, 6.0 ) ;
    rrv_np_base_par -> setVal( 0. ) ;
    rrv_np_base_par -> setConstant( kFALSE ) ;
    allNuisances -> add( *rrv_np_base_par ) ;

    char vname[1000] ;
    sprintf( vname, "prim_mean_%s", NP_base_name ) ;
    RooRealVar* g_mean = new RooRealVar( vname, vname, 0.0,-10.,10. ) ;
    g_mean->setConstant(kTRUE);
    sprintf( vname, "prim_sigma_%s", NP_base_name ) ;
    RooConstVar* g_sigma = new RooConstVar( vname, vname, 1.0 ) ;

    char pdfname[100] ;
    sprintf( pdfname, "pdf_prim_%s", NP_base_name ) ;
    printf("\n\n makeCorrelatedLognormalConstraint : creating base nuisance parameter pdf - %s\n\n", pdfname ) ;
    RooGaussian* base_np_pdf = new RooGaussian( pdfname, pdfname, *rrv_np_base_par, *g_mean, *g_sigma ) ;

    allNuisancePdfs -> add( *base_np_pdf ) ;
    globalObservables -> add( *g_mean ) ;

  }

  //-- create const variables for mean and sigma so that they can be saved and accessed from workspace later.


  char vname[1000] ;
  sprintf( vname, "mean_%s", NP_name ) ;
  RooConstVar* ln_mean  = new RooConstVar( vname, vname, NP_val ) ;

  sprintf( vname, "sigma_%s", NP_name ) ;
  RooConstVar* ln_sigma = new RooConstVar( vname, vname, NP_err ) ;


  RooAbsReal* rar(0x0) ;


  char formula[1000] ;

  if ( !changeSign ) {
    sprintf( formula, "@0 * pow( ( @1/@0 + 1.), @2 )" ) ;
  } else {
    sprintf( formula, "@0 * pow( ( @1/@0 + 1.), -1.0 * @2 )" ) ;
  }

  rar = new RooFormulaVar( NP_name, formula, RooArgSet( *ln_mean, *ln_sigma, *rrv_np_base_par ) ) ;

  printf(" makeCorrelatedLognormalConstraint : creating correlated log-normal NP with formula : %s,  %s, val = %g, mean=%g, sigma=%g\n", formula, NP_name, rar->getVal(), NP_val, NP_err ) ;


  return rar ;

} // makeCorrelatedLognormalConstraint.








//==============================================================================================================

//-- convention for each line in systematics file is
//
//    m1 m2  sig_m1_nt sig_m1_2b sig_m1_3b sig_m1_4b   sig_m2_nt sig_m2_2b sig_m2_3b sig_m2_4b   sig_m3_nt sig_m3_2b sig_m3_3b sig_m3_4b   sig_m4_nt sig_m4_2b sig_m4_3b sig_m4_4b     sb_m1_nt sb_m1_2b sb_m1_3b sb_m1_4b   sb_m2_nt sb_m2_2b sb_m2_3b sb_m2_4b   sb_m3_nt sb_m3_2b sb_m3_3b sb_m3_4b   sb_m4_nt sb_m4_2b sb_m4_3b sb_m4_4b
//
//    where sys is the fractional uncertainty on the signal efficiency for that bin.
//    where the array index, counting from zero, is 2 + sig_or_sb * (Nmet*Nb) + met_index * (Nb) + nb_index
//
//       where sig_or_sb is : =0 for higgs mass signal box, =1 for higgs mass sideband
//             met_index is : =0 for bin1, =1 for bin2, =2 for bin3, =3 for bin4
//             nb_index is  : =0 for no-tag, =1 for 2b, =2 for 3b, =3 for 4b
//
//
bool setupShapeSyst( const char* infile,
                     const char* systname,
                     int constraintType,
                     double target_mgl, double target_mlsp,
                     RooWorkspace& workspace
                     ) {

  printf("\n\n\n setupShapeSyst :  setting up %s systematic.  Input file %s\n\n", systname, infile ) ;

  if ( constraintType == 1 ) {
    printf(" setupShapeSyst : Constraint type for %s : Gaussian (1).\n\n", systname ) ;
  } else if ( constraintType == 2 ) {
    printf(" setupShapeSyst : Constraint type for %s : log-normal (2).\n\n", systname ) ;
  } else {
    printf("  *** setupShapeSyst : Constraint type %d not implemented.\n\n", constraintType ) ;
    return false ;
  }

  char command[1000] ;

  sprintf( command, "head -1 %s | awk '{print NF}'", infile ) ;
  const char* nfields_str = gSystem->GetFromPipe( command ) ;
  int nfields ;
  sscanf( nfields_str, "%d", &nfields ) ;
  printf(" setupShapeSyst: Nfields in %s is %d\n", infile, nfields ) ;


  int ArraySize = nfields ;

  ifstream infq ;
  infq.open(infile) ;
  if ( !infq.good() ) {
    printf("\n\n *** setupShapeSyst: Problem opening input file: %s.\n\n", infile ) ;
    return false ;
  }

  double matchArrayContent[ArraySize] ;
  double nearbyMatchArrayContent[ArraySize] ;

  bool foundMatch = false ;
  bool foundNearbyMatch = false ;


  while ( infq.good() ) {


    double readArrayContent[ArraySize] ;
    for ( int i=0; infq && i<ArraySize; ++ i) {
      infq >> readArrayContent[i] ;
    }

    double mgl  = readArrayContent[0] ;
    double mlsp = readArrayContent[1] ;

    if ( fabs( mgl-target_mgl ) < 1. && fabs( mlsp - target_mlsp ) < 1. ) {

      for ( int i=0; i<ArraySize; i++ ) { matchArrayContent[i] = readArrayContent[i] ; }
      foundMatch = true ;
      printf("\n\n setupShapeSyst : %s :Found mgl=%.0f, mlsp=%.0f\n\n", systname, mgl, mlsp ) ;
      break ;

    }

    if ( fabs( mgl-target_mgl ) < 12.5 && fabs( mlsp - target_mlsp ) < 12.5 && !foundNearbyMatch ) {

      for ( int i=0; i<ArraySize; i++ ) { nearbyMatchArrayContent[i] = readArrayContent[i] ; }
      foundNearbyMatch = true ;
      printf("\n\n setupShapeSyst : %s : Found nearby match mgl=%.0f, mlsp=%.0f  (requested mgl=%.0f, mlsp=%.0f)\n\n", systname, mgl, mlsp, target_mgl, target_mlsp ) ;

    }

  } // reading file?


  double ArrayContent[ArraySize] ;

  if ( foundMatch ) {
    for ( int i=0; i<ArraySize; i++ ) { ArrayContent[i] = matchArrayContent[i] ; }
  } else if ( foundNearbyMatch ) {
    for ( int i=0; i<ArraySize; i++ ) { ArrayContent[i] = nearbyMatchArrayContent[i] ; }
  } else {
    printf("\n\n *** setupShapeSyst : %s : Did not find match or nearby match for mgl=%.0f, mlsp=%.0f\n\n", systname, target_mgl, target_mlsp ) ;
    return false ;
  }


  double syst_msig[3][max_bins_of_met] ;
  double syst_msb[3][max_bins_of_met]  ;

  double minSyst_msig(999.) ;
  double maxSyst_msig(0.) ;
  double minSyst_msb(999.) ;
  double maxSyst_msb(0.) ;

  for ( int smcnbi=0; smcnbi<sigmc_bins_of_nb; smcnbi++ ) {
    for ( int mbi=first_met_bin_array_index; mbi<bins_of_met; mbi++ ) {
      syst_msig[smcnbi][mbi] = ArrayContent[ 2  +  mbi * (sigmc_bins_of_nb)  +  smcnbi ] ;
      syst_msb[smcnbi][mbi]  = ArrayContent[ 2  +  mbi * (sigmc_bins_of_nb)  +  smcnbi  + (bins_of_met*sigmc_bins_of_nb) ] ;
      if ( syst_msig[smcnbi][mbi] > maxSyst_msig ) maxSyst_msig = syst_msig[smcnbi][mbi] ;
      if ( syst_msig[smcnbi][mbi] < minSyst_msig ) minSyst_msig = syst_msig[smcnbi][mbi] ;
      if ( syst_msb[smcnbi][mbi] > maxSyst_msb ) maxSyst_msb = syst_msb[smcnbi][mbi] ;
      if ( syst_msb[smcnbi][mbi] < minSyst_msb ) minSyst_msb = syst_msb[smcnbi][mbi] ;
    } // mbi.
  } // smcnbi.

  printf("\n\n") ;
  printf(" ====== Shape systematics for %s, sig observables\n", systname ) ;
  for ( int smcnbi=0; smcnbi<sigmc_bins_of_nb; smcnbi++ ) {
    printf("  %s, sig, %s :   ", systname, btag_catname[smcnbi] ) ;
    for ( int mbi=first_met_bin_array_index; mbi<bins_of_met; mbi++ ) {
      printf("  %6.3f  ", syst_msig[smcnbi][mbi] ) ;
    } // mbi.
    printf("\n") ;
  } // smcnbi.

  printf("\n\n") ;
  printf(" ====== Shape systematics for %s, sb observables\n", systname ) ;
  for ( int smcnbi=0; smcnbi<sigmc_bins_of_nb; smcnbi++ ) {
    printf("  %s, sb,  %s :   ", systname, btag_catname[smcnbi] ) ;
    for ( int mbi=first_met_bin_array_index; mbi<bins_of_met; mbi++ ) {
      printf("  %6.3f  ", syst_msb[smcnbi][mbi] ) ;
    } // mbi.
    printf("\n") ;
  } // smcnbi.


  printf("\n\n setupShapeSyst: %s : higgs mass sig bins: Min syst = %6.3f, Max syst = %6.3f\n\n", systname, minSyst_msig, maxSyst_msig ) ;
  printf("\n\n setupShapeSyst: %s : higgs mass sb  bins: Min syst = %6.3f, Max syst = %6.3f\n\n", systname, minSyst_msb , maxSyst_msb  ) ;


  for ( int smcnbi=0; smcnbi<sigmc_bins_of_nb; smcnbi++ ) {
    for ( int mbi=first_met_bin_array_index; mbi<bins_of_met; mbi++ ) {

      char pname[100] ;
      bool changeSign ;
      RooAbsReal* rar_par(0x0) ;

      sprintf( pname, "%s_msig_met%d_%s", systname, mbi+1, btag_catname[smcnbi] ) ;
      if ( syst_msig[mbi][smcnbi] < 0 ) { changeSign = true ; } else { changeSign = false ; }
      if ( constraintType == 1 ) {
        rar_par = makeCorrelatedGaussianConstraint(  pname, 1.0, fabs(syst_msig[mbi][smcnbi]), systname, changeSign ) ;
      } else if ( constraintType == 2 ) {
        rar_par = makeCorrelatedLognormalConstraint( pname, 1.0, fabs(syst_msig[mbi][smcnbi]), systname, changeSign ) ;
      } else {
        printf("\n\n *** unknown constraint type : %d\n\n", constraintType ) ;
        return false ;
      }
      cout << flush ;
      workspace.import( *rar_par ) ;

      sprintf( pname, "%s_msb_met%d_%s", systname, mbi+1, btag_catname[smcnbi] ) ;
      if ( syst_msb[mbi][smcnbi] < 0 ) { changeSign = true ; } else { changeSign = false ; }
      if ( constraintType == 1 ) {
        rar_par = makeCorrelatedGaussianConstraint(  pname, 1.0, fabs(syst_msb[mbi][smcnbi]), systname, changeSign ) ;
      } else if ( constraintType == 2 ) {
        rar_par = makeCorrelatedLognormalConstraint( pname, 1.0, fabs(syst_msb[mbi][smcnbi]), systname, changeSign ) ;
      } else {
        printf("\n\n *** unknown constraint type : %d\n\n", constraintType ) ;
        return false ;
      }
      cout << flush ;
      workspace.import( *rar_par ) ;


    } // mbi.
  } // smcnbi.

  printf("\n\n Finished setupShapeSyst for %s.\n\n", systname ) ;
  fflush(stdout) ;

  return true ;

} // setupShapeSyst


//===============================================================================================================

bool readSignalCounts( const char* susy_counts_file, float sig_mass, RooWorkspace& workspace, const char* zh_counts_file, float chi_to_H_bf ) {

  ifstream infp ;
  infp.open( susy_counts_file ) ;
  if ( !infp.good() ) {
    printf("\n\n *** Problem opening input file: %s.\n\n", susy_counts_file ) ;
    return false ;
  }

  int ArraySize = 4 + bins_of_met * sigmc_bins_of_nb * 2 * 2 ; // 2 (sig vs sb) * 2 (vals and errs).

  char command[10000] ;
  sprintf(command, "tail -1 %s | awk '{print NF}' | grep -q %d", susy_counts_file, ArraySize ) ;
  int returnStat = gSystem->Exec(command ) ;
  if ( returnStat !=0 ) {
    printf("\n\n\n *** setSusyScanPoint : expecting %d fields per line in input file %s.  Found ", ArraySize, susy_counts_file ) ; cout << flush ;
    sprintf( command, "tail -1 %s | awk '{print NF}'", susy_counts_file ) ;
    gSystem->Exec(command ) ; cout << flush ;
    printf("\n\n") ;
    return false ;
  }

  float ArrayContent[ArraySize] ;

  bool found(false) ;

  //-- Read in HH signal counts.
  while ( infp.good() ) {

    TString line ;
    line.ReadLine( infp ) ;
    TObjArray* tokens = line.Tokenize(" ") ;
    printf(" number of fields : %d\n", tokens->GetEntries() ) ;
    if ( tokens->GetEntries() != ArraySize ) continue ;

    for ( int i = 0 ; infp && i < ArraySize; ++ i ) {
      TObjString* str = (TObjString*) (tokens->At(i)) ;
      sscanf( (str->GetString()).Data(), "%f", &(ArrayContent[i]) )  ;
    }
    printf( " 1st column: %.0f\n", ArrayContent[0] ) ;

    if ( TMath::Nint( fabs( ArrayContent[0] - sig_mass ) ) == 0   && TMath::Nint( fabs( ArrayContent[1] - 1. ) ) < 2 ) {

      printf("  found mass point %.0f in file %s\n", sig_mass, susy_counts_file ) ;

      found = true ;

      break ;

    }

  } // still reading HH input file?

  if ( !found ) { printf("\n\n *** Did not find mass point %.0f in %s.\n\n", sig_mass, susy_counts_file ) ; return false ; }

  float smc_msig_val[sigmc_bins_of_nb][max_bins_of_met] ;
  float smc_msb_val[sigmc_bins_of_nb][max_bins_of_met] ;

  float smc_msig_err[sigmc_bins_of_nb][max_bins_of_met] ;
  float smc_msb_err[sigmc_bins_of_nb][max_bins_of_met] ;

  for ( int smcnbi=0; smcnbi<sigmc_bins_of_nb; smcnbi++ ) {
    for ( int mbi=first_met_bin_array_index; mbi<bins_of_met; mbi++ ) {

      smc_msig_val[smcnbi][mbi] = (chi_to_H_bf*chi_to_H_bf) * (ArrayContent[ 4 + (                mbi)*(sigmc_bins_of_nb) + smcnbi ]) ;
      smc_msb_val[smcnbi][mbi]  = (chi_to_H_bf*chi_to_H_bf) * (ArrayContent[ 4 + (1*bins_of_met + mbi)*(sigmc_bins_of_nb) + smcnbi ]) ;
      smc_msig_err[smcnbi][mbi] = (chi_to_H_bf*chi_to_H_bf) * (ArrayContent[ 4 + (2*bins_of_met + mbi)*(sigmc_bins_of_nb) + smcnbi ]) ;
      smc_msb_err[smcnbi][mbi]  = (chi_to_H_bf*chi_to_H_bf) * (ArrayContent[ 4 + (3*bins_of_met + mbi)*(sigmc_bins_of_nb) + smcnbi ]) ;

    } // mbi.
  } // smcnbi.

  printf("\n\n\n") ;
  printf("  chi to HH signal counts\n" ) ;
  printf("=============================================================================================================================================================================\n") ;
  printf("  METsig   |        4bSB              4bSIG         |       3bSB              3bSIG         |        2bSB              2bSIG        |        NTSB              NTSIG        |\n") ;
  printf("=============================================================================================================================================================================\n") ;
  fflush(stdout) ;
  for ( int mbi=first_met_bin_array_index; mbi<bins_of_met; mbi++ ) {
    printf(" met bin %d : ", mbi+1 ) ;
    for ( int smcnbi=sigmc_bins_of_nb-1; smcnbi>=0; smcnbi-- ) {
      printf( "  %6.2f +/- %4.2f,  %6.2f +/- %4.2f    |",
              smc_msb_val[smcnbi][mbi] , smc_msb_err[smcnbi][mbi],
              smc_msig_val[smcnbi][mbi], smc_msig_err[smcnbi][mbi] ) ;
    } // smcnbi.
    printf("\n") ;
  } // mbi.
  printf("=====================================================================================================================================\n") ;








  //-- Read in and add ZH counts, if chi_to_H BF is not 1 and have an input file.
  if ( chi_to_H_bf < 1.0 && strlen( zh_counts_file ) > 0 ) {

    ifstream infpzh ;
    infpzh.open( zh_counts_file ) ;
    if ( !infpzh.good() ) {
      printf("\n\n *** Problem opening input file: %s.\n\n", zh_counts_file ) ;
      return false ;
    }

    sprintf(command, "tail -1 %s | awk '{print NF}' | grep -q %d", zh_counts_file, ArraySize ) ;
    int returnStatzh = gSystem->Exec(command ) ;
    if ( returnStatzh !=0 ) {
      printf("\n\n\n *** setSusyScanPoint : expecting %d fields per line in input file %s.  Found ", ArraySize, zh_counts_file ) ; cout << flush ;
      sprintf( command, "tail -1 %s | awk '{print NF}'", zh_counts_file ) ;
      gSystem->Exec(command ) ; cout << flush ;
      printf("\n\n") ;
      return false ;
    }

    found = false ;

    //-- Read in ZH signal counts.
    while ( infpzh.good() ) {

      TString line ;
      line.ReadLine( infpzh ) ;
      TObjArray* tokens = line.Tokenize(" ") ;
      printf(" number of fields : %d\n", tokens->GetEntries() ) ;
      if ( tokens->GetEntries() != ArraySize ) continue ;

      for ( int i = 0 ; infp && i < ArraySize; ++ i ) {
        TObjString* str = (TObjString*) (tokens->At(i)) ;
        sscanf( (str->GetString()).Data(), "%f", &(ArrayContent[i]) )  ;
      }
      printf( " 1st column: %.0f\n", ArrayContent[0] ) ;

      if ( TMath::Nint( fabs( ArrayContent[0] - sig_mass ) ) == 0   && TMath::Nint( fabs( ArrayContent[1] - 1. ) ) < 2 ) {

        printf("  found mass point %.0f in file %s\n", sig_mass, zh_counts_file ) ;

        found = true ;

        break ;

      }

    } // still reading HH input file?

    if ( !found ) { printf("\n\n *** Did not find mass point %.0f in %s.\n\n", sig_mass, zh_counts_file ) ; return false ; }

    for ( int smcnbi=0; smcnbi<sigmc_bins_of_nb; smcnbi++ ) {
      for ( int mbi=first_met_bin_array_index; mbi<bins_of_met; mbi++ ) {

        smc_msig_val[smcnbi][mbi] += 2*(chi_to_H_bf*(1.-chi_to_H_bf)) * (ArrayContent[ 4 + (                mbi)*(sigmc_bins_of_nb) + smcnbi ]) ;
        smc_msb_val[smcnbi][mbi]  += 2*(chi_to_H_bf*(1.-chi_to_H_bf)) * (ArrayContent[ 4 + (1*bins_of_met + mbi)*(sigmc_bins_of_nb) + smcnbi ]) ;
        smc_msig_err[smcnbi][mbi] = sqrt( ::pow( 2*(chi_to_H_bf*(1.-chi_to_H_bf)) * (ArrayContent[ 4 + (2*bins_of_met + mbi)*(sigmc_bins_of_nb) + smcnbi ]), 2.)
                                          + ::pow( smc_msig_err[smcnbi][mbi], 2.) );
        smc_msb_err[smcnbi][mbi]  = sqrt( ::pow( 2*(chi_to_H_bf*(1.-chi_to_H_bf)) * (ArrayContent[ 4 + (3*bins_of_met + mbi)*(sigmc_bins_of_nb) + smcnbi ]), 2.)
                                          + ::pow( smc_msb_err[smcnbi][mbi], 2.) ) ;

      } // mbi.
    } // smcnbi.

    printf("\n\n\n") ;
    printf("  chi to HH+ZH signal counts\n" ) ;
    printf("=============================================================================================================================================================================\n") ;
    printf("  METsig   |        4bSB              4bSIG         |       3bSB              3bSIG         |        2bSB              2bSIG        |        NTSB              NTSIG        |\n") ;
    printf("=============================================================================================================================================================================\n") ;
    fflush(stdout) ;
    for ( int mbi=first_met_bin_array_index; mbi<bins_of_met; mbi++ ) {
      printf(" met bin %d : ", mbi+1 ) ;
      for ( int smcnbi=sigmc_bins_of_nb-1; smcnbi>=0; smcnbi-- ) {
        printf( "  %6.2f +/- %4.2f,  %6.2f +/- %4.2f    |",
                smc_msb_val[smcnbi][mbi] , smc_msb_err[smcnbi][mbi],
                smc_msig_val[smcnbi][mbi], smc_msig_err[smcnbi][mbi] ) ;
      } // smcnbi.
      printf("\n") ;
    } // mbi.
    printf("=====================================================================================================================================\n") ;

  } // reading in and adding ZH counts.









  //-- input stat errors are now in same units as counts (i.e. events).
  //   convert errors into fractional errors.
  for ( int mbi=first_met_bin_array_index; mbi<bins_of_met; mbi++ ) {
    for ( int smcnbi=sigmc_bins_of_nb-1; smcnbi>=0; smcnbi-- ) {
      if ( smc_msig_val[smcnbi][mbi] > 0. ) { smc_msig_err[smcnbi][mbi] = smc_msig_err[smcnbi][mbi] / smc_msig_val[smcnbi][mbi] ; }
      if ( smc_msb_val[smcnbi][mbi]  > 0. ) { smc_msb_err[smcnbi][mbi]  = smc_msb_err[smcnbi][mbi]  / smc_msb_val[smcnbi][mbi]  ; }
    } // smcnbi.
  } // mbi.

  for ( int smcnbi=0; smcnbi<sigmc_bins_of_nb; smcnbi++ ) {
    for ( int mbi=first_met_bin_array_index; mbi<bins_of_met; mbi++ ) {

      char pname[1000] ;

      RooRealVar* rv_smc_msig(0x0) ;
      RooRealVar* rv_smc_msb(0x0) ;
      RooAbsReal* rv_smc_msig_mcstat_syst(0x0) ;
      RooAbsReal* rv_smc_msb_mcstat_syst(0x0) ;

      sprintf( pname, "smc_%s_msig_met%d", btag_catname[smcnbi], mbi+1 ) ;
      rv_smc_msig = new RooRealVar( pname, pname, 0., 1.e6 ) ;
      rv_smc_msig -> setVal( smc_msig_val[smcnbi][mbi] ) ;
      rv_smc_msig -> setConstant( kTRUE ) ;
      workspace.import( *rv_smc_msig ) ;

      sprintf( pname, "smc_%s_msb_met%d", btag_catname[smcnbi], mbi+1 ) ;
      rv_smc_msb = new RooRealVar( pname, pname, 0., 1.e6 ) ;
      rv_smc_msb -> setVal( smc_msb_val[smcnbi][mbi] ) ;
      rv_smc_msb -> setConstant( kTRUE ) ;
      workspace.import( *rv_smc_msb ) ;

      if ( syst_type == 1 ) {
        sprintf( pname, "syst_sig_eff_mc_stats_msig_met%d_%s", mbi+1, btag_catname[smcnbi] ) ;
        rv_smc_msig_mcstat_syst = makeGaussianConstraint( pname, 1.0, smc_msig_err[smcnbi][mbi] ) ;
        sprintf( pname, "syst_sig_eff_mc_stats_msb_met%d_%s", mbi+1, btag_catname[smcnbi] ) ;
        rv_smc_msb_mcstat_syst = makeGaussianConstraint( pname, 1.0, smc_msb_err[smcnbi][mbi] ) ;
      } else if ( syst_type == 2 ) {
        sprintf( pname, "syst_sig_eff_mc_stats_msig_met%d_%s", mbi+1, btag_catname[smcnbi] ) ;
        rv_smc_msig_mcstat_syst = makeLognormalConstraint( pname, 1.0, smc_msig_err[smcnbi][mbi] ) ;
        sprintf( pname, "syst_sig_eff_mc_stats_msb_met%d_%s", mbi+1, btag_catname[smcnbi] ) ;
        rv_smc_msb_mcstat_syst = makeLognormalConstraint( pname, 1.0, smc_msb_err[smcnbi][mbi] ) ;
      } else {
        printf("\n\n *** I don't know how to do syst_type %d\n\n", syst_type ) ; return false ;
      }
      workspace.import( *rv_smc_msig_mcstat_syst ) ;
      workspace.import( *rv_smc_msb_mcstat_syst ) ;

    } // mbi.
  } // smcnbi.

  return true ;


} // readSignalCounts.

//===============================================================================================================

bool setupBtagSFFracMatrix( const char* infracfile, float sig_mass, const char* postfix, RooWorkspace& workspace ) {

  float fracmatrix_p1s[4][4] ;
  float fracmatrix_m1s[4][4] ;

  //-- read in the numbers.
  ifstream instr ;
  instr.open( infracfile ) ;
  if ( !instr.good() ) { printf( "\n\n *** Problem opening %s\n\n", infracfile ) ; return false ; }

  int array_size = 2 + 4*4 + 4*4 ;
  float array_content[ array_size ] ;
  bool found_signal_point(false) ;
  while ( instr.good() ) {
    TString line ;
    line.ReadLine( instr ) ;
    TObjArray* tokens = line.Tokenize( " " ) ;
    if ( tokens -> GetEntries() != array_size ) {
      printf("\n\n *** setupBtagSFFracMatrix : warning.  found %d tokens when expecting %d in %s\n", tokens -> GetEntries(), array_size, infracfile ) ;
      continue ;
    }
    for ( int i = 0 ; instr && i < array_size; ++ i ) {
      TObjString* str = (TObjString*) (tokens->At(i)) ;
      sscanf( (str->GetString()).Data(), "%f", &(array_content[i]) )  ;
    }
    if ( TMath::Nint( fabs( array_content[0] - sig_mass ) ) == 0 ) { found_signal_point = true ; break ; }
  } // still reading?

  if ( ! found_signal_point ) { printf("\n\n *** setupBtagSFFracMatrix : did not find signal point %.0f in %s\n\n", sig_mass, infracfile ) ; return false ; }

  for ( int fci=0; fci<4; fci++ ) {
    for ( int tci=0; tci<4; tci++ ) {
      fracmatrix_p1s[fci][tci] = array_content[ 2 + 4*fci + tci ] ;
      fracmatrix_m1s[fci][tci] = array_content[ 2 + 4*fci + tci + 4*4 ] ;
    } // tci.
  } // fci.





  //-- print them out for the log file.

  printf("\n\n btag SF transfer fraction matrices from %s for signal point %.0f\n\n", infracfile, sig_mass ) ;

  for ( int tci=3; tci>=0; tci-- ) {
    printf("   to %s : ", btag_catname[tci] ) ;
    for ( int fci=0; fci<4; fci++ ) { printf("  %.3f  ", fracmatrix_p1s[fci][tci] ) ; }
    printf("            ") ;
    for ( int fci=0; fci<4; fci++ ) { printf("  %.3f  ", fracmatrix_m1s[fci][tci] ) ; }
    printf("\n") ;
  } // tci.
  printf( "   from    ") ;
  for ( int fci=0; fci<4; fci++ ) { printf("    %s   ", btag_catname[fci] ) ; }
  printf("            ") ;
  for ( int fci=0; fci<4; fci++ ) { printf("    %s   ", btag_catname[fci] ) ; }
  printf("\n\n\n") ;



  //-- check for consistency.
  for ( int fci=0; fci<4; fci++ ) {
    float psum(0.) ;
    float msum(0.) ;
    for ( int tci=0; tci<4; tci++ ) {
      psum += fracmatrix_p1s[fci][tci] ;
      msum += fracmatrix_m1s[fci][tci] ;
    } // tci.
    if ( fabs( 1.0 - psum ) > 0.005 ) { printf(" \n\n *** readBtagSFFracMatrix : fractions for %s, +1 sigma don't add up to 1.\n\n", btag_catname[fci] ) ; return false ; }
    if ( fabs( 1.0 - msum ) > 0.005 ) { printf(" \n\n *** readBtagSFFracMatrix : fractions for %s, -1 sigma don't add up to 1.\n\n", btag_catname[fci] ) ; return false ; }
  } // fci.


  printf("\n\n btag SF transfer fractions successfully read from %s and checked.\n\n", infracfile ) ;



  //-- Look for btag SF nuisance parameter.  If it doesn't exist, make it.

  char prim_name[1000] ;
  sprintf( prim_name, "prim_btag_SF_syst" ) ;
  RooRealVar* rrv_btag_SF_base_par = (RooRealVar*) allNuisances -> find( prim_name ) ;

  if ( rrv_btag_SF_base_par == 0x0 ) {

    printf( " setupBtagSFFracMatrix : creating btag SF nuisance parameter %s\n\n", prim_name ) ;
    rrv_btag_SF_base_par = new RooRealVar( prim_name, prim_name, -6.0, 6.0 ) ;
    rrv_btag_SF_base_par -> setVal( 0. ) ;
    rrv_btag_SF_base_par -> setConstant( kFALSE ) ;
    allNuisances -> add( *rrv_btag_SF_base_par ) ;

    char vname[1000] ;
    sprintf( vname, "mean_%s", prim_name ) ;
    RooRealVar* g_mean = new RooRealVar( vname, vname, 0.0,-10.,10. ) ;
    g_mean->setConstant(kTRUE);
    sprintf( vname, "sigma_%s", prim_name ) ;
    RooConstVar* g_sigma = new RooConstVar( vname, vname, 1.0 ) ;

    char pdfname[100] ;
    sprintf( pdfname, "pdf_%s", prim_name ) ;
    printf("\n\n setupBtagSFFracMatrix : creating btag SF nuisance parameter pdf - %s\n\n", pdfname ) ;
    RooGaussian* base_np_pdf = new RooGaussian( pdfname, pdfname, *rrv_btag_SF_base_par, *g_mean, *g_sigma ) ;

    allNuisancePdfs -> add( *base_np_pdf ) ;
    globalObservables -> add( *g_mean ) ;

  }


  //-- Create the fraction parameters.

  for ( int fci=0; fci<4; fci++ ) {
    for ( int tci=0; tci<4; tci++ ) {

      char pname[100] ;
      RooAbsReal* rar_par ;

      bool same_initial_and_final(false) ;
      if ( fci == tci ) { same_initial_and_final = true ; }

      sprintf( pname, "frac_p1s_from_%s_to_%s_%s", btag_catname[fci], btag_catname[tci], postfix ) ;
      RooRealVar* frac_p1s = new RooRealVar( pname, pname, fracmatrix_p1s[fci][tci], 0., 1. ) ;
      frac_p1s -> setConstant( kTRUE ) ;

      sprintf( pname, "frac_m1s_from_%s_to_%s_%s", btag_catname[fci], btag_catname[tci], postfix ) ;
      RooRealVar* frac_m1s = new RooRealVar( pname, pname, fracmatrix_m1s[fci][tci], 0., 1. ) ;
      frac_m1s -> setConstant( kTRUE ) ;

      //-- RooAsymAbsProd
      //
      //     Returns |g|*f_i,j             if  i != j (same_initial_and_final = false, or off-diagonal matrix element)
      //     Returns (1 - |g|(1 - f_i,j))  if  i == j (same_initial_and_final = true, or on-diagonal matrix element)
      //
      sprintf( pname, "frac_from_%s_to_%s_%s", btag_catname[fci], btag_catname[tci], postfix ) ;
      rar_par = new RooAsymAbsProd( pname, pname, *frac_m1s, *frac_p1s, *rrv_btag_SF_base_par, same_initial_and_final ) ;
      rar_par -> Print("all") ;
      double val_nom, val_p1s, val_m1s ;
      rrv_btag_SF_base_par -> setVal( 0. ) ;
      val_nom = rar_par -> getVal() ;
      rrv_btag_SF_base_par -> setVal( 1. ) ;
      val_p1s = rar_par -> getVal() ;
      rrv_btag_SF_base_par -> setVal(-1. ) ;
      val_m1s = rar_par -> getVal() ;
      printf("    Test for %s (-1s, nom, +1s):  %.3f, %.3f, %.3f\n\n", pname, val_m1s, val_nom, val_p1s ) ;
      rrv_btag_SF_base_par -> setVal(0. ) ;
      workspace.import( *rar_par ) ;

    } // tci.
  } // fci.



  return true ;

} // setupBtagSFFracMatrix

//===============================================================================================================

void makeSystHists( RooWorkspace& workspace, const char* syst_name, const char* syst_var_par, const char* outfile ) {

  printf( "\n\n --- In makeSystHists for %s ------------------------\n\n", syst_name ) ; fflush(stdout) ;

  TFile hfile( outfile, "update" ) ;

  RooRealVar* sig_strength = (RooRealVar*) workspace.obj( "sig_strength" ) ;
  if ( sig_strength == 0x0 ) { printf("\n\n *** makeSystHists : missing sig_strength\n\n" ) ; return ; }
  sig_strength -> setVal( 1.0 ) ;

  for ( int mbi=first_met_bin_array_index; mbi<bins_of_met; mbi++ ) {

    char hname[1000] ;
    char htitle[1000] ;

    sprintf( hname, "h_syst_%s_msig_met%d_nom", syst_name, mbi+1 ) ;
    sprintf( htitle, "Nbtag cat, %s syst, SIG, METsig%d, nominal", syst_name, mbi+1 ) ;
    TH1F* hist_msig_nom = new TH1F( hname, htitle, 4, 0.5, 4.5 ) ;

    sprintf( hname, "h_syst_%s_msig_met%d_m1s", syst_name, mbi+1 ) ;
    sprintf( htitle, "Nbtag cat, %s syst, SIG, METsig%d, -1 sigma", syst_name, mbi+1 ) ;
    TH1F* hist_msig_m1s = new TH1F( hname, htitle, 4, 0.5, 4.5 ) ;

    sprintf( hname, "h_syst_%s_msig_met%d_p1s", syst_name, mbi+1 ) ;
    sprintf( htitle, "Nbtag cat, %s syst, SIG, METsig%d, +1 sigma", syst_name, mbi+1 ) ;
    TH1F* hist_msig_p1s = new TH1F( hname, htitle, 4, 0.5, 4.5 ) ;


    sprintf( hname, "h_syst_%s_msb_met%d_nom", syst_name, mbi+1 ) ;
    sprintf( htitle, "Nbtag cat, %s syst, SB, METsig%d, nominal", syst_name, mbi+1 ) ;
    TH1F* hist_msb_nom = new TH1F( hname, htitle, 4, 0.5, 4.5 ) ;

    sprintf( hname, "h_syst_%s_msb_met%d_m1s", syst_name, mbi+1 ) ;
    sprintf( htitle, "Nbtag cat, %s syst, SB, METsig%d, -1 sigma", syst_name, mbi+1 ) ;
    TH1F* hist_msb_m1s = new TH1F( hname, htitle, 4, 0.5, 4.5 ) ;

    sprintf( hname, "h_syst_%s_msb_met%d_p1s", syst_name, mbi+1 ) ;
    sprintf( htitle, "Nbtag cat, %s syst, SB, METsig%d, +1 sigma", syst_name, mbi+1 ) ;
    TH1F* hist_msb_p1s = new TH1F( hname, htitle, 4, 0.5, 4.5 ) ;


    RooRealVar* syst_par = (RooRealVar*) workspace.obj( syst_var_par ) ;
    if ( syst_par == 0x0 ) { printf("\n\n *** makeSystHists : missing syst par %s\n\n", syst_var_par ) ; return ; }

    for ( int smcnbi=0; smcnbi<sigmc_bins_of_nb; smcnbi++ ) {

      char pname[1000] ;
      double val ;

      sprintf( pname, "mu_sig_%s_msig_met%d", btag_catname[smcnbi], mbi+1 ) ;
      RooFormulaVar* mu_sig = (RooFormulaVar*) workspace.obj( pname ) ;
      printf("\n  %s\n", pname ) ; fflush( stdout ) ;
      if ( mu_sig == 0x0 ) { printf("\n\n *** makeSystHists : missing %s\n\n", pname ) ; return ; }
      mu_sig -> Print() ;

      syst_par -> setVal( 0. ) ;
      val = mu_sig -> getVal() ;
      printf(" nom val %.2f\n", val ) ; fflush( stdout ) ;
      hist_msig_nom -> SetBinContent( smcnbi+1, val ) ;
      hist_msig_nom -> GetXaxis() -> SetBinLabel( smcnbi+1, btag_catname[smcnbi] ) ;

      syst_par -> setVal( -1. ) ;
      val = mu_sig -> getVal() ;
      printf(" m1s val %.2f\n", val ) ; fflush( stdout ) ;
      hist_msig_m1s -> SetBinContent( smcnbi+1, val ) ;
      hist_msig_m1s -> GetXaxis() -> SetBinLabel( smcnbi+1, btag_catname[smcnbi] ) ;

      syst_par -> setVal( +1. ) ;
      val = mu_sig -> getVal() ;
      printf(" p1s val %.2f\n", val ) ; fflush( stdout ) ;
      hist_msig_p1s -> SetBinContent( smcnbi+1, val ) ;
      hist_msig_p1s -> GetXaxis() -> SetBinLabel( smcnbi+1, btag_catname[smcnbi] ) ;



      sprintf( pname, "mu_sig_%s_msb_met%d", btag_catname[smcnbi], mbi+1 ) ;
      RooFormulaVar* mu_sb = (RooFormulaVar*) workspace.obj( pname ) ;
      printf("\n %s\n", pname ) ; fflush( stdout ) ;
      if ( mu_sb == 0x0 ) { printf("\n\n *** makeSystHists : missing %s\n\n", pname ) ; return ; }
      mu_sb -> Print() ;

      syst_par -> setVal( 0. ) ;
      val = mu_sb -> getVal() ;
      printf(" nom val %.2f\n", val ) ; fflush( stdout ) ;
      hist_msb_nom -> SetBinContent( smcnbi+1, val ) ;
      hist_msb_nom -> GetXaxis() -> SetBinLabel( smcnbi+1, btag_catname[smcnbi] ) ;

      syst_par -> setVal( -1. ) ;
      val = mu_sb -> getVal() ;
      printf(" m1s val %.2f\n", val ) ; fflush( stdout ) ;
      hist_msb_m1s -> SetBinContent( smcnbi+1, val ) ;
      hist_msb_m1s -> GetXaxis() -> SetBinLabel( smcnbi+1, btag_catname[smcnbi] ) ;

      syst_par -> setVal( +1. ) ;
      val = mu_sb -> getVal() ;
      printf(" p1s val %.2f\n", val ) ; fflush( stdout ) ;
      hist_msb_p1s -> SetBinContent( smcnbi+1, val ) ;
      hist_msb_p1s -> GetXaxis() -> SetBinLabel( smcnbi+1, btag_catname[smcnbi] ) ;




    } // smcnbi.

    //-- reset it to the nominal val.
    syst_par -> setVal( 0. ) ;

    hist_msig_nom -> Write() ;
    hist_msig_m1s -> Write() ;
    hist_msig_p1s -> Write() ;

    hist_msb_nom -> Write() ;
    hist_msb_m1s -> Write() ;
    hist_msb_p1s -> Write() ;

  } // mbi.

  hfile.Close() ;

} // makeSystHists

//===============================================================================================================

void makeSystHistsMultiPar( RooWorkspace& workspace, const char* syst_name, int npars, const char syst_var_pars[][100], const char* outfile ) {


  printf( "\n\n --- In makeSystHistsMultiPar for %s ------------------------\n\n", syst_name ) ; fflush(stdout) ;

  for ( int pi=0; pi<npars; pi++ ) {
    printf( "   %s : %s\n", syst_name, syst_var_pars[pi] ) ;
  }

  TFile hfile( outfile, "update" ) ;

  RooRealVar* sig_strength = (RooRealVar*) workspace.obj( "sig_strength" ) ;
  if ( sig_strength == 0x0 ) { printf("\n\n *** makeSystHistsMultiPar : missing sig_strength\n\n" ) ; return ; }
  sig_strength -> setVal( 1.0 ) ;

  for ( int mbi=first_met_bin_array_index; mbi<bins_of_met; mbi++ ) {

    char hname[1000] ;
    char htitle[1000] ;

    sprintf( hname, "h_syst_%s_msig_met%d_nom", syst_name, mbi+1 ) ;
    sprintf( htitle, "Nbtag cat, %s syst, SIG, METsig%d, nominal", syst_name, mbi+1 ) ;
    TH1F* hist_msig_nom = new TH1F( hname, htitle, 4, 0.5, 4.5 ) ;

    sprintf( hname, "h_syst_%s_msig_met%d_m1s", syst_name, mbi+1 ) ;
    sprintf( htitle, "Nbtag cat, %s syst, SIG, METsig%d, -1 sigma", syst_name, mbi+1 ) ;
    TH1F* hist_msig_m1s = new TH1F( hname, htitle, 4, 0.5, 4.5 ) ;

    sprintf( hname, "h_syst_%s_msig_met%d_p1s", syst_name, mbi+1 ) ;
    sprintf( htitle, "Nbtag cat, %s syst, SIG, METsig%d, +1 sigma", syst_name, mbi+1 ) ;
    TH1F* hist_msig_p1s = new TH1F( hname, htitle, 4, 0.5, 4.5 ) ;


    sprintf( hname, "h_syst_%s_msb_met%d_nom", syst_name, mbi+1 ) ;
    sprintf( htitle, "Nbtag cat, %s syst, SB, METsig%d, nominal", syst_name, mbi+1 ) ;
    TH1F* hist_msb_nom = new TH1F( hname, htitle, 4, 0.5, 4.5 ) ;

    sprintf( hname, "h_syst_%s_msb_met%d_m1s", syst_name, mbi+1 ) ;
    sprintf( htitle, "Nbtag cat, %s syst, SB, METsig%d, -1 sigma", syst_name, mbi+1 ) ;
    TH1F* hist_msb_m1s = new TH1F( hname, htitle, 4, 0.5, 4.5 ) ;

    sprintf( hname, "h_syst_%s_msb_met%d_p1s", syst_name, mbi+1 ) ;
    sprintf( htitle, "Nbtag cat, %s syst, SB, METsig%d, +1 sigma", syst_name, mbi+1 ) ;
    TH1F* hist_msb_p1s = new TH1F( hname, htitle, 4, 0.5, 4.5 ) ;


    RooRealVar* syst_par[100] ;
    for ( int pi=0; pi<npars; pi++ ) {
      syst_par[pi] = (RooRealVar*) workspace.obj( syst_var_pars[pi] ) ;
      if ( syst_par[pi] == 0x0 ) { printf("\n\n *** makeSystHistsMultiPar : missing syst par %s\n\n", syst_var_pars[pi] ) ; }
    }

    for ( int smcnbi=0; smcnbi<sigmc_bins_of_nb; smcnbi++ ) {

      char pname[1000] ;
      double val ;

      sprintf( pname, "mu_sig_%s_msig_met%d", btag_catname[smcnbi], mbi+1 ) ;
      RooFormulaVar* mu_sig = (RooFormulaVar*) workspace.obj( pname ) ;
      printf("\n  %s\n", pname ) ; fflush( stdout ) ;
      if ( mu_sig == 0x0 ) { printf("\n\n *** makeSystHistsMultiPar : missing %s\n\n", pname ) ; return ; }
      mu_sig -> Print() ;

      for ( int pi=0; pi<npars; pi++ ) { if ( syst_par[pi] != 0 ) { syst_par[pi] -> setVal( 0. ) ; } }
      val = mu_sig -> getVal() ;
      printf(" nom val %.2f\n", val ) ; fflush( stdout ) ;
      hist_msig_nom -> SetBinContent( smcnbi+1, val ) ;
      hist_msig_nom -> GetXaxis() -> SetBinLabel( smcnbi+1, btag_catname[smcnbi] ) ;

      for ( int pi=0; pi<npars; pi++ ) { if ( syst_par[pi] != 0 ) { syst_par[pi] -> setVal( -1. ) ; } }
      val = mu_sig -> getVal() ;
      printf(" m1s val %.2f\n", val ) ; fflush( stdout ) ;
      hist_msig_m1s -> SetBinContent( smcnbi+1, val ) ;
      hist_msig_m1s -> GetXaxis() -> SetBinLabel( smcnbi+1, btag_catname[smcnbi] ) ;

      for ( int pi=0; pi<npars; pi++ ) { if ( syst_par[pi] != 0 ) { syst_par[pi] -> setVal( +1. ) ; } }
      val = mu_sig -> getVal() ;
      printf(" p1s val %.2f\n", val ) ; fflush( stdout ) ;
      hist_msig_p1s -> SetBinContent( smcnbi+1, val ) ;
      hist_msig_p1s -> GetXaxis() -> SetBinLabel( smcnbi+1, btag_catname[smcnbi] ) ;



      sprintf( pname, "mu_sig_%s_msb_met%d", btag_catname[smcnbi], mbi+1 ) ;
      RooFormulaVar* mu_sb = (RooFormulaVar*) workspace.obj( pname ) ;
      printf("\n %s\n", pname ) ; fflush( stdout ) ;
      if ( mu_sb == 0x0 ) { printf("\n\n *** makeSystHistsMultiPar : missing %s\n\n", pname ) ; return ; }
      mu_sb -> Print() ;

      for ( int pi=0; pi<npars; pi++ ) { if ( syst_par[pi] != 0 ) { syst_par[pi] -> setVal( 0. ) ; } }
      val = mu_sb -> getVal() ;
      printf(" nom val %.2f\n", val ) ; fflush( stdout ) ;
      hist_msb_nom -> SetBinContent( smcnbi+1, val ) ;
      hist_msb_nom -> GetXaxis() -> SetBinLabel( smcnbi+1, btag_catname[smcnbi] ) ;

      for ( int pi=0; pi<npars; pi++ ) { if ( syst_par[pi] != 0 ) { syst_par[pi] -> setVal( -1 ) ; } }
      val = mu_sb -> getVal() ;
      printf(" m1s val %.2f\n", val ) ; fflush( stdout ) ;
      hist_msb_m1s -> SetBinContent( smcnbi+1, val ) ;
      hist_msb_m1s -> GetXaxis() -> SetBinLabel( smcnbi+1, btag_catname[smcnbi] ) ;

      for ( int pi=0; pi<npars; pi++ ) { if ( syst_par[pi] != 0 ) { syst_par[pi] -> setVal( +1 ) ; } }
      val = mu_sb -> getVal() ;
      printf(" p1s val %.2f\n", val ) ; fflush( stdout ) ;
      hist_msb_p1s -> SetBinContent( smcnbi+1, val ) ;
      hist_msb_p1s -> GetXaxis() -> SetBinLabel( smcnbi+1, btag_catname[smcnbi] ) ;




    } // smcnbi.

    //-- reset them to their nominal values.
    for ( int pi=0; pi<npars; pi++ ) { if ( syst_par[pi] != 0 ) { syst_par[pi] -> setVal( 0. ) ; } }

    hist_msig_nom -> Write() ;
    hist_msig_m1s -> Write() ;
    hist_msig_p1s -> Write() ;

    hist_msb_nom -> Write() ;
    hist_msb_m1s -> Write() ;
    hist_msb_p1s -> Write() ;

  } // mbi.

  hfile.Close() ;


} // makeSystHistsMultiPar



