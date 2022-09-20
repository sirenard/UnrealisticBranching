#ifndef __SCIP_NODESEL_RANDOM_H__
#define __SCIP_NODESEL_RANDOM_H__

#include <string>
#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#include "objscip/objscip.h"


class NodeSelRandom: public scip::ObjNodesel{
   public:
   explicit NodeSelRandom(
      SCIP* scip
      )
      : ObjNodesel(scip, "random", "Random node selection", 0, 10000)
      {}

   ~NodeSelRandom() override {}


   SCIP_DECL_NODESELSELECT(scip_select) override;
   SCIP_DECL_NODESELCOMP(scip_comp) override;

private:
    static std::string  get_node_hierarchy(SCIP_NODE* node);
}; /*lint !e1712*/



#endif
