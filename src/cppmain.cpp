#include <iostream>
#include <string>

/* include SCIP components */
#include "objscip/objscipdefplugins.h"

#include "Utils.h"
#include "FeaturesCalculator.h"
#include "RegressionModel.h"


using namespace scip;
using namespace std;

/** creates and runs a SCIP instance with default and TSP plugins */
static
SCIP_RETCODE runSCIP(
   int                   argc,               /**< number of arguments from the shell */
   char**                argv                /**< array of shell arguments */
   )
{
   SCIP* scip = nullptr;

    /*********
     * Setup *
     *********/
    Utils::create_scip_instance(&scip, 1);

    SCIP_CALL( SCIPprocessShellArguments(scip, 1, nullptr, "scip.set") );
   /********************
    * Deinitialization *
    ********************/

   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

/** main method starting TSP code */
int main(
   int                   argc,               /**< number of arguments from the shell */
   char**                argv                /**< array of shell arguments */
   )
{
   SCIP_RETCODE retcode;
   srand (time(NULL));

   retcode = runSCIP(argc, argv);
   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
