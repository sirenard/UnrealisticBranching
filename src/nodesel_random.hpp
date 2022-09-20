#ifndef __SCIP_NODESEL_RANDOM_H__
#define __SCIP_NODESEL_RANDOM_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#include "objscip/objscip.h"

class NodeselRandom: public scip::ObjNodesel{
   public:
   NodeselRandom(
      SCIP* scip
      )
      : ObjNodesel(scip, "random", "Random node selection", 0, 10000)
      {}

   virtual ~NodeselRandom() {}

   /*virtual SCIP_DECL_NODESELINIT(nodeselInit)=0;
   virtual SCIP_DECL_NODESELCOPY(nodeselCopyRandom)=0;
   virtual SCIP_DECL_NODESELEXIT(nodeselExit)=0;
   virtual SCIP_DECL_NODESELINITSOL(nodeselInitSol)=0;
   virtual SCIP_DECL_NODESELEXITSOL(nodeselExitSol)=0;*/
   SCIP_DECL_NODESELSELECT(nodeselSelectRandom);
   SCIP_DECL_NODESELCOMP(nodeselCompDfs);
}; /*lint !e1712*/



#endif
