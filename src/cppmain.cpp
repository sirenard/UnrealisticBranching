#include <iostream>
#include <string>

/* include SCIP components */
#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"

#include "nodesel_random.h"
#include "branch_unrealistic.h"
#include "Utils.h"
#include "FeaturesCalculator.h"


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
   assert(argc==2);
   const char* filename = argv[1];

    /*********
     * Setup *
     *********/
    Utils::create_scip_instance(&scip, 1);


    SCIP_CALL(SCIPreadProb(
      scip,
      filename,
      NULL
      )
   );

    /*DatasetWriter writer("node.csv", "../branch.csv");
    FeaturesCalculator featuresCalculator(scip, 1, 1, 1);
    writer.setFeaturesCalculator(&featuresCalculator);
    Branch_unrealistic* branchingRule = (Branch_unrealistic*)SCIPfindBranchrule(scip, "unrealistic");
    Branch_unrealistic::setDataWriter(&writer);*/

    SCIP_CALL( SCIPsolve(scip) );

    SCIP_Longint score = SCIPgetNNodes(scip);
    SCIPinfoMessage(scip, NULL, ("Solved in " + std::to_string(score) + " nodes\n").c_str());


    //SCIP_CALL( SCIPprintBestSol(scip, NULL, FALSE) );

    //SCIP_CALL( SCIPprocessShellArguments(scip, 1, nullptr, "scip.set") );
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
    //srand (time(NULL));
    srand (12);

   retcode = runSCIP(argc, argv);
   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
