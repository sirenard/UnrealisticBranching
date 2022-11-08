#include <iostream>
#include <string>

/* include SCIP components */
#include "objscip/objscipdefplugins.h"

#include "Utils.h"
#include "FeaturesCalculator.h"
#include "RegressionModel.h"
#include "mpi/Worker.h"

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

    SCIP_CALL( SCIPprocessShellArguments(scip, 1, nullptr, "scipmain.set") );
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
    int rank, world_size;
    Worker* node;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    node = new Worker(rank); // if only 1 CPU, the unique node is not seen as a master
    Worker::setInstance(node);

    if(rank == 0){ // master node
        if(world_size > 1)
            node->setWorkersRange(1, world_size);

        SCIP_RETCODE retcode;
        srand (time(NULL));

        retcode = runSCIP(argc, argv);
        if( retcode != SCIP_OKAY )
        {
            SCIPprintError(retcode);
            return -1;
        }
    } else{
        node->work();
    }

    delete node;
    MPI_Finalize();
    return 0;
}
