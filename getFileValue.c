
#include "TSystem.h"
#include <iostream>

//#include "getFileValue.h"

//=====================================================================================================

bool getFileValue( const char* inFile,
                   const char* parameterName,
                   float& returnValue ) {


  returnValue = 1.0 ;

  // Include a blank space at the end to avoid multiple
  //  matches in grep when target is a substring found in other lines.

  char command[10000] ;
  sprintf( command, "grep \"%s \" %s\n", parameterName, inFile ) ;
  TString commandOutput = gSystem->GetFromPipe( command ) ;

  /// printf( " Output of command is : %s\n", commandOutput.Data() ) ;

  char label[1000] ;
  float value ;
  sscanf( commandOutput.Data(), "%s %g", label, &value ) ;
  if ( strcmp( label, parameterName ) == 0 ) {
    printf(" Found %s.  Value is %g\n", parameterName, value ) ;
    returnValue = value ;
    return true ;
  }

  printf("\n\n *** Could not find parameter %s in file %s.\n\n", parameterName, inFile ) ;

  return false ;

}


//=====================================================================================================

bool getFileValueWithError( const char* inFile,
                            const char* parameterName,
                            float& returnValue, float& returnError ) {


  returnValue = 1.0 ;
  returnError = 0.0 ;

  char command[10000] ;
  sprintf( command, "grep \"%s \" %s\n", parameterName, inFile ) ;
  TString commandOutput = gSystem->GetFromPipe( command ) ;

  /// printf( " Output of command is : %s\n", commandOutput.Data() ) ;

  char label[1000] ;
  float value, err ;
  sscanf( commandOutput.Data(), "%s %g %g", label, &value, &err ) ;
  if ( strcmp( label, parameterName ) == 0 ) {
    printf(" Found %s.  Value is %g +/- %g\n", parameterName, value, err ) ;
    returnValue = value ;
    returnError = err ;
    return true ;
  }

  printf("\n\n *** Could not find parameter %s in file %s.\n\n", parameterName, inFile ) ;

  return false ;

}

//=====================================================================================================

bool getFileStringValue( const char* inFile,
                         const char* parameterName,
                         char* returnValue ) {


  sprintf( returnValue, " " ) ;

  // Include a blank space at the end to avoid multiple
  //  matches in grep when target is a substring found in other lines.

  char command[10000] ;
  sprintf( command, "grep \"%s \" %s\n", parameterName, inFile ) ;
  TString commandOutput = gSystem->GetFromPipe( command ) ;

  /// printf( " Output of command is : %s\n", commandOutput.Data() ) ;

  char label[1000] ;
  char the_string[10000] ;
  sscanf( commandOutput.Data(), "%s %s", label, the_string ) ;
  if ( strcmp( label, parameterName ) == 0 ) {
    printf(" Found %s.  Value is %s\n", parameterName, the_string ) ;
    sprintf( returnValue, "%s", the_string ) ;
    return true ;
  }

  printf("\n\n *** Could not find parameter %s in file %s.\n\n", parameterName, inFile ) ;

  return false ;

}


//=====================================================================================================


bool getFileMultiStringValue( const char* inFile,
                              const char* parameterName,
                              int& nvals,
                              char returnVals[][1000] ) {



  // Include a blank space at the end to avoid multiple
  //  matches in grep when target is a substring found in other lines.

  char command[10000] ;
  sprintf( command, "grep \"%s \" %s\n", parameterName, inFile ) ;
  TString commandOutput = gSystem->GetFromPipe( command ) ;

  /// printf( " Output of command is : %s\n", commandOutput.Data() ) ;

  char label[1000] ;
  sscanf( commandOutput.Data(), "%s", label ) ;
  if ( strcmp( label, parameterName ) == 0 ) {
    printf(" Found %s.  \n", parameterName ) ;

    TObjArray* strings = commandOutput.Tokenize(" ") ;
    printf(" Breaks into %d tokens.\n", strings -> GetEntries() ) ;
    nvals = strings -> GetEntries() - 1 ;
    for ( int i=0; i<nvals; i++ ) {
      TObjString* str = (TObjString*) (strings->At(i+1)) ;
      printf( " string %d : %s\n", i, str->GetString().Data() ) ;
      sprintf( returnVals[i], "%s", str->GetString().Data() ) ;
    }

    return true ;
  }

  printf("\n\n *** Could not find parameter %s in file %s.\n\n", parameterName, inFile ) ;

  return false ;

}


//=====================================================================================================












