#define SCIP_DEBUG
#include "nodesel_random.h"
#include "scip/pub_message.h"
#include "scip/pub_nodesel.h"
#include "scip/pub_tree.h"
#include "scip/scip_message.h"
#include "scip/scip_tree.h"
#include <cstring>
#include <stdlib.h>
#include <iostream>

#define NODESEL_NAME             "random"
#define NODESEL_DESC             "random selection"
#define NODESEL_STDPRIORITY           0
#define NODESEL_MEMSAVEPRIORITY  100000



/*
 * Callback methods
 */

/** node selection method of node selector */
SCIP_DECL_NODESELSELECT(NodeSelRandom::scip_select)
{  /*lint --e{715}*/
   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);
   assert(selnode != NULL);


    *selnode = SCIPgetBestNode(scip);
    //SCIPcutoffNode(scip, *selnode);

    //SCIPdebugMsg(scip, ("Selection: " + std::to_string(    SCIPnodeGetNumber(*selnode)) + "\n").c_str());
    //std::string msg = get_node_hierarchy(*selnode) + "\n";
    //SCIPdebugMsg(scip, msg.c_str());
    //SCIPprintNodeRootPath(scip, *selnode, NULL);

    //SCIPnodeIsPropagatedAgain(*selnode);



    /*if(selnode == NULL){
        *selnode = SCIP
    }*/

   return SCIP_OKAY;
}


/** node comparison method of node selector */
SCIP_DECL_NODESELCOMP(NodeSelRandom::scip_comp)
{
   int res = rand()%2-1;
    //return SCIPnodeGetNumber(node1) < SCIPnodeGetNumber(node2)?1:-1;
    //SCIPdebugMsg(scip, (std::to_string(res) + "\n").c_str());
   return res;
}

std::string NodeSelRandom::get_node_hierarchy(SCIP_NODE *node) {
    SCIP_NODE* parent = SCIPnodeGetParent(node);
    int num = SCIPnodeGetNumber(node);

    if (parent==nullptr){
        return std::to_string(num);
    } else{
        return std::to_string(num) + " - " + get_node_hierarchy(parent);
    }
}
